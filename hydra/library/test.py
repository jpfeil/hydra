import cyrand
import numpy as np

X = np.array([1, 1, 1, 2, 2, 2, 3, 3, 3], dtype=np.int64) - 1
Y = np.array([1, 1, 2, 2, 2, 2, 3, 3, 4], dtype=np.int64) - 1

print( cyrand.ri(X, Y) )
