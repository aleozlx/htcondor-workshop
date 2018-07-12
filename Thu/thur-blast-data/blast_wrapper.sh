#!/bin/bash

tar xzf pdbaa_files.tar.gz
./blastx -db pdbaa -query $1 -out $1.result
rm pdbaa.*
