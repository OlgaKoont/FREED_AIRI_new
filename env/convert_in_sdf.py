import os
from openbabel import pybel


input_directory = '/mnt/tank/scratch/okonovalova/freedpp/freedpp/env/docking_results_mor/'
output_directory = input_directory  


os.makedirs(output_directory, exist_ok=True)

for filename in os.listdir(input_directory):
    if filename.endswith('.pdbqt'):
        pdbqt_file = os.path.join(input_directory, filename)
        sdf_file = os.path.join(output_directory, filename.replace('.pdbqt', '.sdf'))
        
        try:
            
            mol = pybel.readfile("pdbqt", pdbqt_file).__next__()  
            
            
            mol.write("sdf", sdf_file, overwrite=True)
            print(f"Converted: {filename} to {os.path.basename(sdf_file)}")
        
        except Exception as e:
            print(f"Error converting {filename}: {e}")
    
    