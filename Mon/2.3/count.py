#!/usr/bin/env python

import os
import sys
import operator

if len(sys.argv) != 2:
    print 'Usage: %s DATA' % (os.path.basename(sys.argv[0]))
    sys.exit(1)
input_filename = sys.argv[1]

words = {}

my_file = open(input_filename, 'r')
for line in my_file:
    line_words = line.split()
    for word in line_words:
        if word in words:
            words[word] += 1
        else:
            words[word] = 1
my_file.close()

sorted_words = sorted(words.items(), key=operator.itemgetter(1))
for word in sorted_words:
    print '%s %8d' % (word[0], word[1])
