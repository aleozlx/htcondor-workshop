#!/bin/bash
wget -S http://stash.osgconnect.net/~username/pdbaa_files.tar.gz
tar xzf pdbaa_files.tar.gz
./blastx -db pdbaa -query $1 -out $1.result
rm pdbaa.*
