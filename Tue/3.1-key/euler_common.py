#!/bin/env python

"""
Anagramic Squares (Euler problem 98): https://projecteuler.net/problem=98
Common functions 
"""
import sys

def usage(args, err_msg=''):
    '''Print error message and usage instructions'''
    if err_msg:
        err_msg = 'ERROR: %s\n' %err_msg
    msg = "%sUsage: %s %s" % (err_msg, sys.argv[0], args)
    sys.exit(msg)

def find_anagrams(candidates):
    '''
    Take list of ints or words and return a list, 
    containing lists of anagrams
    '''
    candidates = [str(x) for x in candidates]
    anagrams = dict()

    # combine all words/ints into groups of anagrams    
    for x in candidates:
        key = ''.join(sorted(list(x)))
        if key not in anagrams.keys():
            anagrams[key] = set([x])
        anagrams[key].add(x)
        
    # drop words/ints without any anagrams
    for k in anagrams.keys():
        if len(anagrams[k]) == 1:
            anagrams.pop(k)

    return anagrams.values()
