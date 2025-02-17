#!/bin/sh
export CUDA_VISIBLE_DEVICES=4
python ../ffreed/main.py \
	    --exp_root experiments \
	    --alert_collections ../../alert_collections.csv \
	    --fragments ../../zinc_crem.json \
            --starting_smile "c1([*:1])c([*:2])ccc([*:3])c1" \
	    --fragmentation crem \
            --device 'cuda' \
            --batch_size 32 \
            --db_device 2 \
	    --n_conf 1 \
            --num_mols 500 \
	    --save_freq 50 \
	    --epochs 200 \
	    --commands "evaluate" \
	    --reward_version soft \
            --seed 42 \
            --weights "1.0" \
            --objectives "DBAffinity" \
            --name db_4djh_42_500_200e \
            --db_protein_pdb ../../DynamicBind/data/4djh_cleaned.pdb \
            --db_header main_db_42_500_200e \
            --db_results result_folder \
            --db_env /nfs/home/okonovalova/anaconda3/envs/dynamicbind/bin/python \
            --db_relax_env /nfs/home/okonovalova/anaconda3/envs/relax/bin/python \
            --checkpoint experiments/db_4djh_42_500_200e/ckpt/model_200.pth
            

           #--checkpoint experiments/db_4djh_42_500_100e/ckpt/model_200.pth
           # --ligandFile ../DynamicBind/ligand.csv \


