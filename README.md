
# Regularized Dynamic Time Warping Kernel (KDTW)

KDTW is a similarity measure dedicated to (multivariate) time serie matching constructed from DTW with the property that KDTW(.,.) is a positive definite kernel (homogeneous to an inner product in the so-called Reproducing Kernel Hilbert Space). Following earlier work by Cuturi & al.  [1], namely the so-called Global Alignment kernel (GA-kernel or GAK), the derivation of KDTW is detailed in Marteau & Gibet 2014  [4]. KDTW is a R-Convolution kernel as defined in  [2]. 

## Algorithmic complexity
The algorithmic complexity of KDTW is O(N^2), where N is the length of the longest of the two time series given in entry, and when no corridor size is specified. This complexity drop down to O(C.N) when a symmetric corridor of size C is exploited. 

## Interpretation
Please see [4] for more details.


## References

[1]   M. Cuturi, J.-P. Vert, Ø. Birkenes, and Tomoko Matsui. A kernel for time series based on global alignments. In Acoustics, Speech and Signal Processing, 2007. ICASSP 2007. IEEE International Conference on, volume 2, pages II–413–II–416, April 2007.

[2]   David Haussler. Convolution kernels on discrete structures. Technical Report UCS-CRL-99-10, University of California at Santa Cruz, Santa Cruz, CA, USA, 1999.

[3]   E. J. Keogh, X. Xi, L. Wei, and C.A. Ratanamahatana. The UCR time series classification-clustering datasets, 2006. http://wwwcs.ucr.edu/ eamonn/time_series_data/.

[4]   Pierre-François Marteau and Sylvie Gibet. On Recursive Edit Distance Kernels with Application to Time Series Classification. IEEE Transactions on Neural Networks and Learning Systems, pages 1–14, June 2014. 14 pages.

[5]   Bernhard Schölkopf, Alexander Smola, and Klaus-Robert Müller. Nonlinear component analysis as a kernel eigenvalue problem. Neural Comput., 10(5):1299–1319, July 1998.
