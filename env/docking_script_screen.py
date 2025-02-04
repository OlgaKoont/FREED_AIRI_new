import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import meeko
from meeko import PDBQTWriterLegacy

smiles_file = 'kor_mor_jamda_info.csv'
receptor_file = 'protein_4djh_chain_A.pdbqt'
output_directory = 'docking_results_kor/'
vina_program = '/mnt/tank/scratch/okonovalova/freedpp/freedpp/env/qvina02'


os.makedirs(output_directory, exist_ok=True)

df = pd.read_csv(smiles_file)
smiles_list = df['Smiles'].tolist()  


vina_config = 'vina_4djh.cfg'

for i, smiles in enumerate(smiles_list):
    ligand = Chem.MolFromSmiles(smiles)
    protonated_ligand = Chem.AddHs(ligand)
    Chem.AllChem.EmbedMolecule(protonated_ligand)

    molecule_preparation = meeko.MoleculePreparation()
    # Assuming you have already prepared your molecule
    molecule_setups = molecule_preparation.prepare(protonated_ligand)  # This returns a list of MoleculeSetup instances
#    lig_pdbqt = molecule_preparation.write_pdbqt_string()

    writer = PDBQTWriterLegacy()

    # Write each setup to a PDBQT file
    for setup in molecule_setups:
        pdbqt_string = writer.write_string(setup)
        with open("ligand.pdbqt", "w") as f:
            f.write(str(pdbqt_string))
   
    ligand_pdbqt = os.path.join(output_directory, f'ligand_{i}.pdbqt')
    output_file = os.path.join(output_directory, f'output_ligand_{i}.pdbqt')
    log_file = os.path.join(output_directory, f'log_ligand_{i}.txt')
    
    # Run docking command
    command = f'{vina_program} --config {vina_config} --receptor {receptor_file} --ligand {ligand_pdbqt} --out {output_file} --log {log_file}'
    os.system(command)