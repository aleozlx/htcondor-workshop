#!/bin/sh

/bin/tar -cvzf perm_combined_$1.tar.gz perm_part_$1_*.Rdat > tar_$1.log
