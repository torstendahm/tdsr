import math
import numpy as np

arg =  [-1000, -100, -30, -10 , -1, 0.001, 1.,  10., 30.,  100., 1000., 10000.]
arg = np.asarray( arg )
narg = len(arg)
print(' a', arg)

print(' Waterlevel fuer failuretime')
argmax = 30.
argnull = 0.
X = np.exp(arg)
print(' X', X)
X = np.exp(arg, out = np.ones(narg), where= (arg < argmax) )
print(' mit where')
print(' X', X)

print(' ')
print(' Waterlevel fuer probability')
X = np.exp(-arg, out = np.ones(narg), where= (-arg < argnull) )
print(' X', X)
