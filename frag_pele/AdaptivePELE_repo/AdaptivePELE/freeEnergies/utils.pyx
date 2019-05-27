import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport sqrt

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] sum_rows(double[:,:] C):
    cdef Py_ssize_t sx, sy
    cdef double tmp_sum
    sx = C.shape[0]
    sy = C.shape[1]
    cdef double[:] x = np.empty((sy))
    for i in range(sx):
        tmp_sum = 0.0
        for j in range(sx):
            tmp_sum += C[i, j]
        x[i] = tmp_sum

    return x


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def buildRevTransitionMatrix(double[:,:] C):
    """
        #Implemented as Prinz paper
    """
    cdef double[:,:] X, T
    cdef Py_ssize_t sx, sy
    sx = C.shape[0]
    sy = C.shape[1]
    X = np.empty((sx, sy))
    T = np.empty((sx, sy))
    for k in range(sx):
        for ky in range(sy):
            X[k, ky]=C[k,ky]+C[ky,k]
            T[k, ky]= 0
    cdef double[:] x, c
    x = sum_rows(X)
    c = sum_rows(C)
    cdef int iterations = 1000
    cdef int it, i, j
    cdef double a, b, z, tmp, tmp2, tmp3

    for it in range(iterations):
        # if it != 0 and it % 10 == 0: print it
        for i in range(sx):
            X[i,i] = C[i,i] * (x[i] - X[i,i]) / (c[i] - C[i,i])
        x = sum_rows(X)
        for i in range(sx - 1):
            for j in range(i + 1, sx):
                a = c[i] - C[i,j] + c[j] - C[j,i]
                b = c[i]*(x[j]-X[i,j]) + c[j]*(x[i] - X[i,j]) -(C[i,j] + C[j,i])*(x[i] + x[j] - 2*X[i,j])
                z = -(C[i,j] + C[j,i])*(x[i] - X[i,j])*(x[j] - X[i,j])
                tmp = b**2 - 4*a*z
                tmp2 = 2*a
                tmp3 = sqrt(tmp)
                X[j,i] = (-b+tmp3)/tmp2
                X[i,j] = X[j,i]
        x = sum_rows(X)

    for i in range(sx):
        for j in range(sx):
            T[i,j] = X[i,j]/x[i]
    return np.array(T)
