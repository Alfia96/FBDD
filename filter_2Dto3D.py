
from rdkit import Chem
from rdkit.Chem import Descriptors

def count_molecules(input_sdf):
    supplier = Chem.SDMolSupplier(input_sdf)
    count = 0
    for mol in supplier:
        if mol is not None:
            count += 1
    return count

def filter_molecules(input_sdf, output_sdf):
    # Count the number of molecules in the input SDF file
    num_molecules = count_molecules(input_sdf)
    print(f"Total number of molecules in the SDF file before filtering based on Rule of Three criteria: {num_molecules}")

    # Read the input SDF file
    supplier = Chem.SDMolSupplier(input_sdf)
    writer = Chem.SDWriter(output_sdf)
    
    count_filtered = 0
    
    for mol in supplier:
        if mol is None:
            continue
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        clogp = Descriptors.MolLogP(mol)
        
        # Apply Rule of Three criteria
        if mw <= 300 and hbd <= 3 and hba <= 3 and clogp <= 3:
            writer.write(mol)
            count_filtered += 1
    
    writer.close()
    
    print(f"Number of molecules that meet the Rule of Three criteria: {count_filtered}")

#Path to input and output file
input_sdf = 'scaffolds_5500.sdf'  
output_sdf = 'Trial_filter.sdf'  

filter_molecules(input_sdf, output_sdf)

