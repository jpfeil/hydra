import numpy as np
cimport numpy as np

from cython.parallel import prange

cdef extern from "../dev/rand_index.h":
    double c_randi(int *x,
                   int maxx,
                   int *y,
                   int maxy,
                   int size);

cpdef binomial(long n, long k):
    cdef long[:] C = np.zeros(k+1, dtype=np.int64)
    C[0] = 1        # nC0 is 1
    cdef long i, j
    for i in range(1, n + 1):
        #print("i: ", i)
        for j in range(min(i, k), 0, -1):
            #print("j: ", j)
            C[j] = C[j] + C[j-1]
    return C[k]

cdef randi(long[:] x, long[:] y):
    cdef long n = x.shape[0]
    cdef int nx = max(x) + 1
    cdef int ny = max(y) + 1
    cdef long[:,:] mat = np.zeros((nx, ny), dtype=np.int64)
    cdef int i, j
    cdef double S1=0.
    cdef double S2=0.
    cdef double S3=0.

    for i in range(x.shape[0]):
        mat[x[i], y[i]] += 1

    ai = np.sum(mat, axis=1)
    bj = np.sum(mat, axis=0)

    for i in range(nx):
        for j in range(ny):
            S1 += binomial(mat[i, j], 2)

    for i in range(nx):
        S2 += binomial(ai[i], 2)

    for j in range(ny):
        S3 += binomial(bj[j], 2)

    bit = (S2 * S3) / binomial(n, 2)
    cdef double ri = (S1 - bit) / (0.5 * (S2 + S3) - bit)

    if np.isclose(ri, 0., rtol=0.1, atol=0.1):
        #print("Close to zero")
        ri = 0.

    assert 0. <= ri <= 1, 'Value is out of bounds! %f' % ri
    return ri