#!/bin/bash

for i in {1..7500}

do
	(
        python ~/Sweetie_hunter/pdb2pqr/pdb2pqr.py --ff=amber $i.pdb $i.tmp #generates charge and coords using python code
	grep -Ev 'MOL' $i.tmp |grep -E 'ATOM|HETATM' | awk '{print $6, $7, $8, $9}' > $i.pqr
	) 
done 


