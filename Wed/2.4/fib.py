#!/usr/bin/env python3
import sys
import os

if len(sys.argv) != 2:
    print('Usage: %s MAXIMUM' % (os.path.basename(sys.argv[0])))
    sys.exit(1)
maximum = int(sys.argv[1])
n1 = n2 = 1
while n2 <= maximum:
    n1, n2 = n2, n1 + n2
print('The greatest Fibonacci number up to %d is %d' % (maximum, n1))

