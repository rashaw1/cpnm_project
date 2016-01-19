import numpy as np

N = 100000
rtN = np.sqrt(N)/2.0
f = open('huge.inp', 'w')
f.write('basis,\nC, 0.2\nbasisend\ngeom,\n')

randoms = np.random.uniform(-rtN, rtN, 3*N)
i = 0
while (i < N):
    f.write('C, %f, %f, %f\n' % (randoms[i], randoms[i+1], randoms[i+2]))
    i += 3

f.write('geomend\nthreshold, 1e-6')
f.close()
