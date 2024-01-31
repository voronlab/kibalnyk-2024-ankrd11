#!/bin/bash
set -euxo pipefail

# create data directory
mkdir -p 'data/mtx'

# fetch counts
curl 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE258835&format=file' -o 'data/mtx/GSE258835.tar'
tar -xvf 'data/mtx/GSE258835.tar' -C 'data/mtx'
rm 'data/mtx/GSE258835.tar'

for f in 'data/mtx/'*'.tar.gz'; do
	exdir="data/mtx/$(basename "$(sed 's/GSM[0-9]*_//' <<<"${f//-mat/}")" '.tar.gz')"
	mkdir "$exdir"
	tar -xvzf "$f" -C "$exdir"
	rm "$f"
done
