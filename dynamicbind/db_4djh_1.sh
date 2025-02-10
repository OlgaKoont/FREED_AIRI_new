#!/bin/sh
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5
python ../main_db_min.py \
	    --exp_root experiments \
	    --alert_collections ../alert_collections.csv \
	    --fragments ../zinc_crem.json \
            --starting_smile "c1([*:1])c([*:2])ccc([*:3])c1" \
	    --fragmentation crem \
            --device 'cuda' \
            --device_num 2 \
	    --n_conf 1 \
            --num_mols 1 \
	    --save_freq 1 \
	    --epochs 1 \
            --savings_per_complex 1 \
            --samples_per_complex 1 \
	    --commands "train,sample" \
	    --reward_version soft \
            --seed 42 \
            --hts \
            --weights "-1.0" \
            --objectives "DB" \
            --name db_4djh_42_32_32_rew \
            --proteinFile ../DynamicBind/data/4djh_cleaned.pdb \
            --savings_per_complex 1 \
            --inference_steps 20 \
            --header main_db_42_32_32_rew \
            --results result_folder \
            --python /nfs/home/okonovalova/anaconda3/envs/dynamicbind/bin/python \
            --relax_python /nfs/home/okonovalova/anaconda3/envs/relax/bin/python 
            

           #--checkpoint experiments/db_4djh_42/ckpt/model_200.pth
           # --ligandFile ../DynamicBind/ligand.csv \


