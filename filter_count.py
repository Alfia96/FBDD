from rdkit import Chem
from rdkit.Chem import Descriptors

def filter_molecules(input_sdf, output_sdf):
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

# Input and output file path
input_sdf = 'scaffolds.sdf'  
output_sdf = 'filtered_molecules_count.sdf'  

filter_molecules(input_sdf, output_sdf)
