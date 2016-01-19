import numpy as np

FILENAME = 'thc.xyz'
NEWFILE = 'thc.inp'

inf = open(FILENAME, 'r')
outf = open(NEWFILE, 'w')
outf.write('geom,\n')

for line in inf:
    line.strip()
    cols = line.split()
    outf.write('%s, %s, %s, %s\n' % (cols[0], cols[1], cols[2], cols[3]))

outf.write('geomend\n')
outf.close()
inf.close()
