cdef extern from "rand_index.h":
    double c_randi(int *x,
                   int maxx,
                   int *y,
                   int maxy,
                   int size);

def crandi(int[:] x, int maxx, int[:] y, int maxy, int size):
    return c_randi(&x[0], maxx, &y[0], maxy, size)
       

