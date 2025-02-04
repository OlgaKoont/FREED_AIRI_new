import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel


smiles_file = 'kor_mor_jamda_info.csv'
receptor_file = 'protein_8efo_chain_R.pdbqt'
output_directory = 'docking_results_mor/'
vina_program = '/mnt/tank/scratch/okonovalova/freedpp/freedpp/env/qvina02'


os.makedirs(output_directory, exist_ok=True)

df = pd.read_csv(smiles_file)
smiles_list = df['Smiles'].tolist()  


vina_config = 'vina_8efo.cfg'

# Assuming smiles_list and output_directory are defined
for i, smiles in enumerate(smiles_list):
    # Create a molecule from SMILES
    ligand = pybel.readstring("smi", smiles)
    
    # Add hydrogens and generate 3D coordinates
    ligand.OBMol.AddHydrogens()
    ligand.make3D()  # Generate 3D coordinates

    # Write the ligand to a temporary PDBQT file
    temp_pdbqt_file = f"temp_ligand_{i}.pdbqt"
    ligand.write("pdbqt", temp_pdbqt_file, overwrite=True)

    # Move the temporary file to the desired output directory with a proper name
    ligand_pdbqt = os.path.join(output_directory, f'ligand_{i}.pdbqt')
    os.rename(temp_pdbqt_file, ligand_pdbqt)

    output_file = os.path.join(output_directory, f'output_ligand_{i}.pdbqt')
    log_file = os.path.join(output_directory, f'log_ligand_{i}.txt')
    
    # Run docking command
    command = f'{vina_program} --config {vina_config} --receptor {receptor_file} --ligand {ligand_pdbqt} --out {output_file} --log {log_file}'
    os.system(command)
    
    # Convert the output PDBQT file to PDB format using Open Babel
    output_pdb_file = os.path.join(output_directory, f'output_ligand_{i}.pdb')
    
    # Read the output PDBQT file and write it as a PDB file
    ob_output_ligand = pybel.readfile("pdbqt", output_file)
    
    for mol in ob_output_ligand:
        mol.write("pdb", output_pdb_file, overwrite=True)

    print("Docking and conversion completed.")