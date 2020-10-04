
# Regularized Dynamic Time Warping Kernel (KDTW)

**kdtw.py** and **kdtw_cdist** are two python3.* implementations of KDTW, a similarity measure dedicated to (multivariate) time serie matching. KDTW is derived from DTW while ensuring the property that KDTW(.,.) is a positive definite kernel (homogeneous to an inner product in the so-called Reproducing Kernel Hilbert Space). Following earlier work by Cuturi & al.  [1], namely the so-called Global Alignment kernel (GA-kernel or GAK), the derivation of KDTW is detailed in Marteau & Gibet 2014 [3]. KDTW is a R-Convolution kernel as defined in [2]. 

## Algorithmic complexity
The algorithmic complexity of KDTW is O(N^2), where N is the length of the longest of the two time series given in entry, and when no corridor size is specified. This complexity drop down to O(C.N) when a symmetric corridor of size C is exploited. 

## Implementations
* kdtw_cdist.py uses numpy and scipy.spatial.distance.cdist to accelerate the computation of the local kernel matrix
* kdtw.py uses only numpy and is slower than kdtw_cdist

## Interpretation
The main idea leading to the construction of positive definite kernels from a given elastic distance used to match two time series is simply the following: instead of keeping only one of the best alignment paths, the new kernel sums up the costs of all the existing sub-sequence alignment paths with some smoothing factor that will favor good alignments while penalizing bad alignments. In addition, this smoothing factor can be optimized. Thus, basically we replace the min (or max) operator by a summation operator and introduce a symmetric corridor function (the function h in the recursive equation above) to optionally limit the summation and the computational cost. Then we add a regularizing recursive term (Kxx) such that the proof of the positive definiteness property can be understood as a direct consequence of the Haussler’s convolution theorem  [2]. Please see  [3] for the details.


## References

[1]   M. Cuturi, J.-P. Vert, Ø. Birkenes, and Tomoko Matsui. A kernel for time series based on global alignments. In Acoustics, Speech and Signal Processing, 2007. ICASSP 2007. IEEE International Conference on, volume 2, pages II–413–II–416, April 2007.

[2]   David Haussler. Convolution kernels on discrete structures. Technical Report UCS-CRL-99-10, University of California at Santa Cruz, Santa Cruz, CA, USA, 1999.

[3]   Pierre-François Marteau and Sylvie Gibet. On Recursive Edit Distance Kernels with Application to Time Series Classification. IEEE Transactions on Neural Networks and Learning Systems, pages 1–14, June 2014. 14 pages.
