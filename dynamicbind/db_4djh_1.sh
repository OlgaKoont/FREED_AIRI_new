#!/bin/sh
python ../main_db.py \
	    --exp_root experiments \
	    --alert_collections ../alert_collections.csv \
	    --fragments ../zinc_crem.json \
            --starting_smile "c1([*:1])c([*:2])ccc([*:3])c1" \
	    --fragmentation crem \
	    --n_conf 1 \
            --num_mols 1 \
	    --save_freq 50 \
	    --epochs 1 \
            --savings_per_complex 1 \
	    --commands "train,sample" \
	    --reward_version soft \
            --seed 42 \
            --weights "1.0" \
            --objectives "DB" \
            --name db_4djh__42_1 \
            --proteinFile ../DynamicBind/data/4djh_cleaned.pdb \
            --savings_per_complex 1 \
            --inference_steps 20 \
            --header main_db_42_1 \
            --results result_folder \
            --python /nfs/home/okonovalova/anaconda3/envs/dynamicbind/bin/python \
            --relax_python /nfs/home/okonovalova/anaconda3/envs/relax/bin/python 
            

           #--checkpoint experiments/db_4djh_42/ckpt/model_200.pth
           # --ligandFile ../DynamicBind/ligand.csv \


