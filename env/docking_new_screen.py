#!/usr/bin/env python
import os
from multiprocessing import Pool
from subprocess import run
import glob
import numpy as np
import argparse
import os
import json
from rdkit import Chem
from utils_first import lmap

def str2floats(s):
    return lmap(float, s.split(','))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Docking
    parser.add_argument('--receptor', required=True)
    parser.add_argument('--box_center', required=True, type=str2floats)
    parser.add_argument('--box_size', required=True, type=str2floats)
    parser.add_argument('--ligand', required=True)
    parser.add_argument('--out')
    parser.add_argument('--vina_program', required=True)
    parser.add_argument('--exhaustiveness', type=int, default=8)
    parser.add_argument('--num_modes', type=int, default=10)
    parser.add_argument('--num_sub_proc', type=int, default=1)
    parser.add_argument('--n_conf', type=int, default=3)
    parser.add_argument('--error_val', type=float, default=99.9)
    parser.add_argument('--timeout_gen3d', type=int, default=None)
    parser.add_argument('--timeout_dock', type=int, default=None)
    parser.add_argument('--seed', help='RNG seed', type=int, default=42)

    return parser.parse_args()


output_dir = "complex_exp"
os.makedirs(output_dir, exist_ok=True)

class DockingVina:
    def __init__(self, config):
        self.config = config

    def __call__(self, smile):
        affinities = []
        for i in range(self.config['n_conf']):
            os.environ['OB_RANDOM_SEED'] = str(self.config['seed'] + i)
            affinities.append(DockingVina.docking(smile, **self.config))
        return min(affinities)

    @staticmethod
    def docking(smile, *, vina_program, receptor, box_center,
                box_size, error_val, seed, num_modes, exhaustiveness,
                timeout_dock, timeout_gen3d, **kwargs):

        
        ligand_file = os.path.join(output_dir, f'ligand_{seed}.pdbqt')
        docking_file = os.path.join(output_dir, f'docking_{seed}.pdbqt')

        
        run_line = f"obabel -:{smile} --gen3D -h -opdbqt -O {ligand_file}"
        result = run(run_line.split(), capture_output=True, text=True,
                     timeout=timeout_gen3d, env=os.environ)

        if "Open Babel Error" in result.stdout or "3D coordinate generation failed" in result.stdout:
            return error_val

        
        run_line = f"{vina_program} --receptor {receptor} --ligand {ligand_file} --out {docking_file}"
        run_line += f" --center_x {box_center[0]} --center_y {box_center[1]} --center_z {box_center[2]}"
        run_line += f" --size_x {box_size[0]} --size_y {box_size[1]} --size_z {box_size[2]}"
        run_line += f" --num_modes {num_modes} --exhaustiveness {exhaustiveness} --seed {seed}"

        result = run(run_line.split(), capture_output=True, text=True, timeout=timeout_dock)

        return DockingVina.parse_output(result.stdout, error_val)

    @staticmethod
    def parse_output(result, error_val):
        result_lines = result.split('\n')
        check_result = False
        affinity = error_val

        for result_line in result_lines:
            if result_line.startswith('-----+'):
                check_result = True
                continue
            if not check_result:
                continue
            if result_line.startswith('Writing output'):
                break
            if result_line.startswith('Refine time'):
                break
            lis = result_line.strip().split()
            if not lis[0].isdigit():
                break
            affinity = float(lis[1])
            break
        return affinity

def main(args):
    DockingVina(args)

def setup():
    args = parse_args()
    return args


if __name__ == '__main__':
    args = setup()
    main(args)
