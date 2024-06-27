import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

#Function to filter molecules based on NRB and TPSA
def filter_nrb_tpsa(molecule):
    if molecule is None:
        return False
    nrb = rdMolDescriptors.CalcNumRotatableBonds(molecule)
    tpsa = rdMolDescriptors.CalcTPSA(molecule)
    return nrb <= 3 and tpsa <= 60

#Paths to input and output SDF files
input_sdf = 'Chembridge_filter_2.sdf'
output_sdf = 'Chembridge_filter_2_final.sdf'

#Initialize the count of molecules
initial_count = 0
filtered_count = 0

#Open the input and output sdf files
supplier = Chem.SDMolSupplier(input_sdf)
writer = Chem.SDWriter(output_sdf)

for mol in supplier:
    initial_count += 1
    if mol is None:
        print(f"Warning: Skipping molecule {initial_count} (invalid structure)")
        continue
    try:
        if filter_nrb_tpsa(mol):
            writer.write(mol)
            filtered_count += 1
    except Exception as e:
        print(f"Error processing molecule {initial_count}: {e}")

#Close the writer
writer.close()

#Print the results
print(f"Total molecules in input file: {initial_count}")
print(f"Total molecules in output file: {filtered_count}")