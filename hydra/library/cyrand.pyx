import numpy as np
from cython.parallel import prange

cdef extern from "rand_index.h":
    double crandi(long *x,
                  int maxx,
                  long *y,
                  int maxy,
                  int size)

def ri(long[:] x, long[:] y):
    cdef int size = np.int32( x.shape[0] )
    cdef int maxx = np.int32( max(x) + 1 )
    cdef int maxy = np.int32( max(y) + 1 )

    return crandi(&x[0], maxx, &y[0], maxy, size)

# Won't compile because of an issue with
# slicing. I will get back to this soon.
#def pri(long[:,:] x, long[:,:] y):
#    cdef int size = np.int32( x.shape[0] )
#    cdef double[:] res = np.zeros(size)
#    cdef int i
#    cdef int[:] maxxs = np.zeros(size)
#    cdef int[:] maxys = np.zeros(size)
#
#    for i in range(size):
#        maxxs[i] = max(x[i, :]) + 1
#        maxys[i] = max(y[i, :]) + 1
#
#
#    for i in prange(size, nogil=True):
#        res[i] = crandi(&x[i][0],
#                        maxxs[i],
#                        &y[i][0],
#                        maxys[i],
#                        size)
#    return res
       

