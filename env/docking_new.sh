#!/bin/bash

python3 docking_new_screen.py \
	    --receptor ../protein_4djh_chain_A.pdbqt \
            --ligand 'C1=CC2c3ccccc3NC(C(Nc3ccc4c(-c5cccc(-c6ccnc7ccccc67)c5)c[nH]c4c3)c3c[nH]c4ccccc34)C2C1' \
            --vina_program /mnt/tank/scratch/okonovalova/freedpp/freedpp/env/qvina02 \
            --exhaustiveness 1 \
            --box_center "4.139,-22.228,60.007" \
    	    --box_size "17.206,25.0,20.574" \
	    --seed 42 \
            --out kor_lig_1

