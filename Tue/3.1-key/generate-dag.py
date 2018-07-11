#!/bin/env python

import sys

input_file = open('find-anagrams.out', 'r')
anagrams = input_file.readlines()
input_file.close()

if len(anagrams) == 0:
    sys.exit('ERROR: find-anagrams.py did not produce any anagrams')

dag_contents = ''
i = 0
for line in anagrams:
    serial = str(i).zfill(3)
    node_name = 'Squares%s' % serial
    anagram_list = line.strip().split(',')
    if len(anagram_list) < 2:
        sys.exit('ERROR: Got args list with fewer than two anagrams')
    anagram_args = ' '.join(anagram_list)
    dag_contents += 'JOB %s squares.sub\n' % node_name
    dag_contents += 'VARS %s anagrams="%s" serial="%s"\n' \
                    % (node_name, anagram_args, serial)
    i += 1

dag_file = open ('squares.dag', 'w')
dag_file.write(dag_contents)
dag_file.close()
