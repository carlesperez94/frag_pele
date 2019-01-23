cimport numpy as np
import numpy as np

def buildRevTransitionMatrix_fast(C, iterations=5):
    """
        #Implemented as Prinz paper
    """
    X = C + C.T
    x = X.sum(axis=1)
    c = C.sum(axis=1)

    return loop(iterations, X, x, C, c)

def loop(int iterations, np.ndarray[double, ndim=2] X, np.ndarray[double, ndim=1] x, np.ndarray[double, ndim=2] C, np.ndarray[double, ndim=1] c):

    cdef int n = len(x)

    for it in range(iterations):
        for i in range(n):
            X[i,i] = C[i,i] * (x[i] - X[i,i]) / (c[i] - C[i,i])

        x = X.sum(axis=1)

        c_vec = c[:,np.newaxis]
        x_vec = x[:,np.newaxis]

        A = c_vec - C + c_vec.T - C.T
        B = c_vec*(x_vec.T - X) + c_vec.T*(x_vec - X) - (C + C.T) * (x_vec + x_vec.T - 2 * X)
        Z = -(C+C.T) * (x_vec - X) * (x_vec.T - X)

        indicesU = np.triu_indices(n,1)
        X[indicesU] = (-B[indicesU] + np.sqrt(B[indicesU]**2 - 4 * A[indicesU]*Z[indicesU]))/ (2 * A[indicesU])
        X.T[indicesU] = X[indicesU]


        x = X.sum(axis=1)

    return X / x[:,np.newaxis]
