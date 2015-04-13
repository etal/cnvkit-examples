#!/bin/bash
for f in $@
do
	mv -v ${f} ${f:5:7}.cw-pair.cns
done
