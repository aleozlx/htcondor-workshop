#!/bin/bash

module load xrootd
module load stashcp
stashcp /user/s18_alexyan/public/$1 ./
./ffmpeg -i $1 -b:v 400k -s 640x360 $2
rm $1

