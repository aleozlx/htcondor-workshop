#!/bin/sh
echo 'Date: ' `date` 
echo 'Host: ' `hostname` 
echo 'System: ' `uname -spo` 
echo "Program: $0" 
echo "Args: $*"
echo 'ls: ' `ls`
# END
