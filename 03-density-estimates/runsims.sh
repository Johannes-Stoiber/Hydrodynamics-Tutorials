#!/bin/bash

n=16
while [[ $n -le 256 ]] ; do
	echo "nbg = $n"
	cat box1.param | sed "s|XXXX|${n}|g" > box.param
	mpirun -np 1 OpenGadget3/P-Gadget3 box.param
	python3 calc_rho.py $n
	(( n = n+8))
done
