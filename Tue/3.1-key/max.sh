#!/bin/bash

# Collate values of anagrammic squares and find their maximum value
cat squares/square-*.out | sort -n | tail -n 1
