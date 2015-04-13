#!/bin/bash
for g in $@
do
	mv -v ${g} ${g:5:7}.cw-noref.cns
done
