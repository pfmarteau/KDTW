import numpy as np
import KDTW

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    A=np.array([[0],[0],[1],[1],[2],[3],[5],[2],[0],[1],[-0.1]], dtype=float)
    B=np.array([[0],[1],[2],[2.5],[3],[3.5],[4],[4.5],[5.5],[2],[0],[0],[.25],[.05],[0]], dtype=float)
    C=np.array([[4],[4],[3],[3],[3],[3],[2],[5],[2],[.5],[.5],[.5]], dtype=float)

    sigma=1.0
    epsilon=1e-20
    print("kdtw(A,B)=", KDTW.similarity(A,B,sigma, epsilon))
    print("kdtw(A,C)=", KDTW.similarity(A,C,sigma, epsilon))
    print("kdtw(B,C)=", KDTW.similarity(B,C,sigma, epsilon))

    plt.plot(A, label='A')
    plt.plot(B, label='B')
    plt.plot(C, label='C')
    plt.legend()
    plt.show()

