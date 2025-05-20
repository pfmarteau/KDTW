#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 16:46:45 2020
@author: pfm
"""
import numpy as np

# Filename: kdtw.py
# Python source code for the "Kernelized" Dynamic Time Warping similarity (as defined in the reference below).
# Author: Pierre-Francois Marteau
# Version: V1.0 du 13/09/2014,
# Licence: GPL
# ******************************************************************
# This software and description is free delivered "AS IS" with no
# guaranties for work at all. Its up to you testing it modify it as
# you like, but no help could be expected from me due to lag of time
# at the moment. I will answer short relevant questions and help as
# my time allow it. I have tested it played with it and found no
# problems in stability or malfunctions so far.
# Have fun.
# *****************************************************************
# Please cite as:
# @article{marteau:hal-00486916,
#   AUTHOR = {Marteau, Pierre-Francois and Gibet, Sylvie},
#   TITLE = {{On Recursive Edit Distance Kernels with Application to Time Series Classification}},
#   JOURNAL = {{IEEE Transactions on Neural Networks and Learning Systems}},
#   PAGES = {1-14},
#   YEAR = {2014},
#   MONTH = Jun,
#   KEYWORDS = {Elastic distance, Time warp kernel, Time warp inner product, Definiteness, Time series classification, SVM},
#   DOI = {10.1109/TNNLS.2014.2333876},
#   URL = {http://hal.inria.fr/hal-00486916}
# }


def kdtw(A, B, sigma=1, epsilon=1e-3):
    """Dynamic programming implementation of KDTW kernel

    :param A: first multivariate time series: array of array (n x d), n is the number
              of samples, d is the dimension of each sample
    :param B: second multivariate time series: array of array (n x d), n is the number
              of samples, d is the dimension of each sample
    :param sigma: must be >0; used in the exponential local kernel
    :param epsilon: must be 1 > epsilon > 0, used in the exponential local kernel
    :return: similarity between A and B (the higher, the more similar)
    """
    n = np.shape(A)[0] + 1
    m = np.shape(B)[0] + 1

    DP = np.zeros((n, m))
    DP1 = np.zeros((n, m))
    DP2 = np.zeros(max(n, m))

    DP2[0] = 1.0
    for i in range(1, min(n, m)):
        DP2[i] = Dlpr(A[i-1], B[i-1], sigma, epsilon)

    DP[0, 0] = 1
    DP1[0, 0] = 1

    for i in range(1, n):
        DP[i, 0] = DP[i - 1, 0] * Dlpr(A[i-1], B[0], sigma, epsilon)
        DP1[i, 0] = DP1[i - 1, 0] * DP2[i]

    for j in range(1, m):
        DP[0, j] = DP[0, j - 1] * Dlpr(A[0], B[j-1], sigma, epsilon)
        DP1[0, j] = DP1[0, j - 1] * DP2[j]

    for i in range(1, n):
        for j in range(1, m):
            lcost = Dlpr(A[i-1], B[j-1], sigma, epsilon)
            DP[i, j] = (DP[i - 1, j] + DP[i, j - 1] + DP[i - 1, j - 1]) * lcost
            DP1[i, j] = DP1[i - 1, j] * DP2[i] + DP1[i, j - 1] * DP2[j]
            if i == j:
                DP1[i, j] += DP1[i - 1, j - 1] * lcost
    DP = DP + DP1
    return DP[n - 1, m - 1]


def Dlpr(a, b, sigma=1, epsilon=1e-3):
    # local similarity between two samples
    # a: 1d numpy array
    # b: 1d numpy array
    # input sigma: >0 used in the exponential local kernel
    # input epsilon: 1 > epsilon > 0
    # return the local matching similarity (probability)
    return (np.exp(-np.sum((a - b) ** 2) / sigma) + epsilon) / (3 * (1 + epsilon))


def compare_to_matlab():
    # parameters
    sigma = 0.125
    epsilon = 1e-20

    # Univariate case
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=float).reshape(-1, 1)
    y = np.array([5, 6, 7, 8, 9, 1, 2], dtype=float).reshape(-1, 1)
    expected_distance = 1.2814254453822292e-102

    distance = kdtw(x, y, sigma=sigma, epsilon=epsilon)
    np.testing.assert_almost_equal(distance, expected_distance, decimal=112)

    # multivariate case
    x = np.array(
        [[1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 2, 1, 2, 3, 2, 1]], dtype=float
    ).T
    y = np.array([[5, 6, 7, 8, 9, 1, 2], [5, 6, 7, 8, 7, 6, 5]], dtype=float).T
    expected_distance = 1.828989485754932e-183

    distance = kdtw(x, y, sigma=sigma, epsilon=epsilon)
    np.testing.assert_almost_equal(distance, expected_distance, decimal=193)


# Simple test
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    compare_to_matlab()

    A = np.array([[0], [0], [1], [1], [2], [3], [5], [2], [0], [1], [-0.1]])
    B = np.array(
        [
            [0],
            [1],
            [2],
            [2.5],
            [3],
            [3.5],
            [4],
            [4.5],
            [5.5],
            [2],
            [0],
            [0],
            [0.25],
            [0.05],
            [0],
        ]
    )
    C = np.array([[4], [4], [3], [3], [3], [3], [2], [5], [2], [0.5], [0.5], [0.5]])

    print("kdtw(A,B)=", kdtw(A, B, 1))
    print("kdtw(A,C)=", kdtw(A, C, 1))
    print("kdtw(B,C)=", kdtw(B, C, 1))

    plt.plot(A, label="A")
    plt.plot(B, label="B")
    plt.plot(C, label="C")
    plt.legend()
    plt.show()
