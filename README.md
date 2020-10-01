<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>
# Regularized Dynamic Time Warping Kernel (KDTW)

KDTW is a similarity measure dedicated to (multivariate) time serie matching constructed from DTW with the property that KDTW(.,.) is a positive definite kernel (homogeneous to an inner product in the so-called Reproducing Kernel Hilbert Space). Following earlier work by Cuturi & al.  [1], namely the so-called Global Alignment kernel (GA-kernel), the derivation of KDTW is detailed in Marteau & Gibet 2014  [4]. KDTW is a convolution kernel as defined in  [2]. After giving a recursive definition of KDTW.

## Definition
The recursive definition of $$K_{DTW}$$ kernel as defined in ~\cite{Marteau2014} is as follows:
\\
\begin{align}
\label{Eq.MEREDK}
\begin{array}{ll}
\mathcal{K}_{DTW}(X_{p},Y_{q})=K^{xy}_{DTW}(X_{p}, Y_{q})+K^{xx}_{DTW}(X_{p},Y_{q}) \\
\\ \text{where} 
\\
K^{xy}_{DTW}(X_{p},Y_{q}) = \beta \cdot e^{- d_{E}^{2}(x(p),y(q))/\sigma}  \\
   \sum \left\{
   \begin{array}{ll}
    h(p-1,q)K^{xy}_{DTW}(X_{p-1},Y_{q}) \\
   h(p-1,q-1) K^{xy}_{DTW}(X_{p-1},Y_{q-1})  \\
    h(p,q-1)K^{xy}_{DTW}(X_{p},Y_{q-1}) \\
   \end{array}
   \right.\\
\\ \text{and} 
\\
   K^{xx}_{DTW}(X_{p},Y_{q}) = \beta \cdot \\
   \sum \left\{
   \begin{array}{ll}
    h(p-1,q) K^{xx}_{DTW}(X_{p-1},Y_{q})e^{- d_{E}^{2}(x(p),y(p))/\sigma}  \\
    \Delta_{p,q} h(p,q)K^{xx}_{DTW}(X_{p-1},Y_{q-1})e^{- d_{E}^{2}(x(p),y(q))/\sigma}   \\
    h(p,q-1)K^{xx}_{DTW}(X_{p},Y_{q-1})e^{-d_{E}^{2}(x(q),y(q))/\sigma} \\
   \end{array}
   \right.\\
  \end{array}
\end{align}
where 

\begin{itemize}
\item $\Delta_{p,q}$ is the Kronecker's symbol, 
\item $h$ is a symmetric binary non negative function, usually in $\{0,1\}$, used to define a symmetric corridor around the main diagonal to limit the "time elasticity" of the kernel,  
\item $\sigma \in \mathbb{R}^{+}$ is a \textit{bandwidth} parameter which weights the local contributions, i.e. the distances between locally aligned positions, 
\item $\beta \in [1/3, 1]$ is a normalization factor. Given an order of magnitude for time series lengths, $\beta$ and $\sigma$ are dependent. Namely if $\sigma$ tends to $\infty$ (is large comparatively to the variance of the time series samples), then $\beta$ will have to tend to $\frac{1}{3}$. If $\sigma$ tends to $0$ (is small comparatively to the variance of the time series samples) then $\beta$ will have to tend to $1$. By default, $\beta=\frac{1}{3}$, which prevents the "infinity explosion" of the kernel.
\item $d_E(.,.)$ is the Euclidean distance defined on $\mathbb{R}^{k}$ (or any negative conditionally definite kernel define on $\mathbb{R}^{k}$).

\end{itemize} 

\section{Algorithmic complexity}
The recursive definition given above shows that the algorithmic complexity of $K_{DTW}$ is $O(N^2)$, where $N$ is the length of the longest of the two time series given in entry, and when no corridor size is specified. This complexity drop down to $O(C.N)$ when a symmetric corridor of size $C$ is exploited. 

\section{Interpretation}
The main idea leading to the construction of positive definite
kernels from a given elastic distance used to match two time series is simply the following: instead of keeping only one of the best alignment paths, the
new kernel will try to sum up the costs of all the existing sub-sequence
alignment paths with some weighting factor that will
favor good alignments while penalizing bad alignments. In
addition, this weighting factor can be optimized.  Thus, basically we replace the \textit{min} (or \textit{max}) operator by a summation operator and introduce a symmetric corridor function (the function $h$ in the recursive equation above) to optionally limit the summation and the computational cost. Then we add a regularizing recursive term ($K^{xx}$) such that the proof of the positive definiteness property can be understood as a direct
consequence of the Haussler’s convolution theorem ~\cite{Haussler99}. Please see ~\cite{Marteau2014} for more details.


References
[1]   M. Cuturi, J.-P. Vert, Ø. Birkenes, and Tomoko Matsui. A kernel for time series based on global alignments. In Acoustics, Speech and Signal Processing, 2007. ICASSP 2007. IEEE International Conference on, volume 2, pages II–413–II–416, April 2007.

[2]   David Haussler. Convolution kernels on discrete structures. Technical Report UCS-CRL-99-10, University of California at Santa Cruz, Santa Cruz, CA, USA, 1999.

[3]   E. J. Keogh, X. Xi, L. Wei, and C.A. Ratanamahatana. The UCR time series classification-clustering datasets, 2006. http://wwwcs.ucr.edu/ eamonn/time_series_data/.

[4]   Pierre-François Marteau and Sylvie Gibet. On Recursive Edit Distance Kernels with Application to Time Series Classification. IEEE Transactions on Neural Networks and Learning Systems, pages 1–14, June 2014. 14 pages.

[5]   Bernhard Schölkopf, Alexander Smola, and Klaus-Robert Müller. Nonlinear component analysis as a kernel eigenvalue problem. Neural Comput., 10(5):1299–1319, July 1998.
