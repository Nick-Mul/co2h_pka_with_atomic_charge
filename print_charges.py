from rdkit import Chem
from rdkit.Chem import AllChem
import sys

def read_mol2_with_charges(file_path):
    # Read the MOL2 file
    mol = Chem.MolFromMol2File(file_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read MOL2 file: {file_path}")
    
    # Ensure the molecule has 3D coordinates
    if not mol.GetConformer().Is3D():
        raise ValueError("Molecule does not have 3D coordinates")
    
    return mol

def find_carboxylic_acid(mol):
    # SMARTS pattern for carboxylic acid
    patt = Chem.MolFromSmarts('C(=O)[OH0]')
    matches = mol.GetSubstructMatches(patt)
    
    if not matches:
        #print("No carboxylic acid group found.")
        return None
    
    # Return the first match (in case there are multiple carboxylic groups)
    return matches[0]

def extract_carboxylic_charges(mol, carboxylic_indices):
    charges = {}
    for idx in carboxylic_indices:
        atom = mol.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        charge = atom.GetDoubleProp('_TriposPartialCharge')  # Assumes Gasteiger charges are present
        charges[f"{symbol}{idx}"] = charge
    return charges

def process_mol2_file(file_path):
    atom_list = []
    file_list = []
    charge_list = []
    try:
        mol = read_mol2_with_charges(file_path)
        carboxylic_indices = find_carboxylic_acid(mol)
        
        if carboxylic_indices:
            charges = extract_carboxylic_charges(mol, carboxylic_indices)
            #print("Partial charges for carboxylic acid group:")
            for atom, charge in charges.items():
                print(f"{atom},{charge:.4f},{file_path}")
                atom_list.append(atom)
                charge_list.append(charge)
                file_list.append(file_path)
        else:
            #print(f"No carboxylic acid group found in the molecule {file_path}.")
            pass
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

    return atom_list, charge_list, file_list

# Usage
file_path = sys.argv[1]
(process_mol2_file(file_path))