from util import *
from nzmath import ffpack_interface
import time
import numpy as np

field = GF(101)
data = []

try:
    for i in range(1, 300, 3):
        print i, "by", i, "matrices"
        ffpack_times = []
        nzmath_times = []
        for j in range(5):
            M1 = random_matrix(i, field)
            M2 = random_matrix(i, field)
            
            ffpack_interface.ENABLE_FFPACK = True
            t = time.time()
            M3 = M1 * M2
            ffpack_times.append(time.time() - t)
            
            ffpack_interface.ENABLE_FFPACK = False
            t = time.time()
            M3 = M1 * M2
            nzmath_times.append(time.time() - t)
        data.append((i, np.mean(ffpack_times), np.std(ffpack_times),
                        np.mean(nzmath_times), np.std(nzmath_times)))
except KeyboardInterrupt:
    pass

f = open("multiply.dat", "w")
print >>f, "n ffpack_mean ffpack_std nzmath_mean nzmath_std"
for entry in data:
    print >>f, "%d %f %f %f %f" % entry
f.close()
