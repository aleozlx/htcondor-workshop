#!/bin/env python

"""
Anagramic Squares (Euler problem 98): https://projecteuler.net/problem=98

This script takes a list of space-separated words as arguments and finds
the largest anagramic number, if any
"""
import euler_common
import io
import random
import sys
import time

if random.randint(0,1):
    f = io.BytesIO()
    f.write(' '*10**7)
    time.sleep(60)

usage_args = '<list of anagrams>'
if len(sys.argv) < 3:
    euler_common.usage('<list of anagrams>',
                       'Anagram list must contain at least two words')

words = [x.lower() for x in sys.argv[1:]]
words_len = set([len(x) for x in words])

if len(words_len) != 1:
    euler_common.usage('<list of anagrams>',
                       'Words must be the same length')
num_digits = words_len.pop()

def set_len(iterable):
    '''Return the number of unique elements in an iterable'''
    return len(set(iterable))

# construct lits of potential squares
x = 1
sq = x**2
squares = []
sq_digits = len(str(sq))
while sq_digits <= num_digits:
    if sq_digits == num_digits:
        squares.append(sq)
    x += 1
    sq = x**2
    sq_digits = len(str(sq))

ana_squares = set()
idx_word = list(words.pop())
for sq in squares:
    str_sq = list(str(sq))
    if set_len(idx_word) != set_len(str_sq):
        continue # Letter->number mapping must be unique
    idx = dict(zip(idx_word, str_sq))
    for wrd in words:
        letters = []
        for char in list(wrd):
            letters.append(idx[char])
        indexed_num = int(''.join(letters))

        if indexed_num in squares:
            ana_squares = ana_squares.union(set([sq, indexed_num]))

try:
    print max(ana_squares)
except ValueError:
    pass
