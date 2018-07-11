#!/bin/env python

"""
Anagramic Squares (Euler problem 98): https://projecteuler.net/problem=98

This script takes a file containing a list of comma-separated words and prints
lists of anagrams to stdout
"""

import euler_common
import sys

args_usage = '<anagrams filename>'
try:
    fname = sys.argv[1]
    f = open(fname, 'r')
    text = f.read()
    f.close()
except IndexError:
    euler_common.usage(args_usage, 'Missing filename argument.')
except IOError, exc:
    euler_common.usage(args_usage, exc)

words = [x.replace('"', '') for x in text.split(',')]        

# combine all words into groups of anagrams
anagrams = dict()
for w in words:
    key = ''.join(sorted(list(w)))
    if key not in anagrams.keys():
        anagrams[key] = [w]
    elif w not in anagrams[key]:
        anagrams[key].append(w)

# drop words without any anagrams
for k in anagrams.keys():
    if len(anagrams[k]) == 1:
        anagrams.pop(k)

for x in anagrams.values():
    print ','.join(x)
