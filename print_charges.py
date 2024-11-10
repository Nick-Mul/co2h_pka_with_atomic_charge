from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os
from typing import Tuple, List, Dict, Optional
import glob
import pandas as pd

def read_mol2_with_charges(file_path: str) -> Chem.Mol:
    """
    Read a MOL2 file and return an RDKit molecule object with charges.

    Args:
        file_path (str): Path to the MOL2 file.

    Returns:
        Chem.Mol: RDKit molecule object.

    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file can't be read as a MOL2 or doesn't have 3D coordinates.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    mol = Chem.MolFromMol2File(file_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read MOL2 file: {file_path}")
    
    if not mol.GetConformer().Is3D():
        raise ValueError("Molecule does not have 3D coordinates")
    
    return mol

def find_carboxylic_acid(mol: Chem.Mol) -> Optional[Tuple[int, ...]]:
    """
    Find the first carboxylic acid group in the molecule.

    Args:
        mol (Chem.Mol): RDKit molecule object.

    Returns:
        Optional[Tuple[int, ...]]: Indices of atoms in the carboxylic acid group, or None if not found.
    """
    patt = Chem.MolFromSmarts('C(=O)[OH0]')
    matches = mol.GetSubstructMatches(patt)
    return matches[0] if matches else None

def extract_carboxylic_charges(mol: Chem.Mol, carboxylic_indices: Tuple[int, ...]) -> Dict[str, float]:
    """
    Extract partial charges for atoms in the carboxylic acid group.

    Args:
        mol (Chem.Mol): RDKit molecule object.
        carboxylic_indices (Tuple[int, ...]): Indices of atoms in the carboxylic acid group.

    Returns:
        Dict[str, float]: Dictionary of atom symbols and their charges.
    """
    return {
        f"{mol.GetAtomWithIdx(idx).GetSymbol()}{idx}": mol.GetAtomWithIdx(idx).GetDoubleProp('_TriposPartialCharge')
        for idx in carboxylic_indices
    }

def process_mol2_file(file_path: str) -> Tuple[List[str], List[float], List[str]]:
    """
    Process a MOL2 file and extract carboxylic acid charges.

    Args:
        file_path (str): Path to the MOL2 file.

    Returns:
        Tuple[List[str], List[float], List[str]]: Lists of atom symbols, charges, and file paths.
    """
    atom_list, charge_list, file_list = [], [], []

    try:
        mol = read_mol2_with_charges(file_path)
        carboxylic_indices = find_carboxylic_acid(mol)
        
        if carboxylic_indices:
            charges = extract_carboxylic_charges(mol, carboxylic_indices)
            for atom, charge in charges.items():
                print(f"{atom},{charge:.4f},{file_path}")
                atom_list.append(atom)
                charge_list.append(charge)
                file_list.append(file_path)
        
    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}", file=sys.stderr)

    return atom_list, charge_list, file_list

def process_directory(directory: str) -> Tuple[List[str], List[float], List[str]]:
    """
    Process all MOL2 files in a given directory.

    Args:
        directory (str): Path to the directory containing MOL2 files.

    Returns:
        Tuple[List[str], List[float], List[str]]: Combined lists of atom symbols, charges, and file paths.
    """
    all_atoms, all_charges, all_files = [], [], []

    # Use glob to get all .mol2 files in the directory
    mol2_files = glob.glob(os.path.join(directory, "*.mol2"))

    for file_path in mol2_files:
        atoms, charges, files = process_mol2_file(file_path)
        all_atoms.extend(atoms)
        all_charges.extend(charges)
        all_files.extend(files)

     # Create DataFrame
    df = pd.DataFrame({
        'Atom': all_atoms,
        'Charge': all_charges,
        'File': all_files
    })

    # Pivot the table
    pivoted_df = df.pivot(index='File', columns='Atom', values='Charge')

    # Sort the columns
    pivoted_df = pivoted_df.reindex(sorted(pivoted_df.columns), axis=1)

    # Reset the index to make 'Molecule' a column
    pivoted_df = pivoted_df.reset_index()

    return pivoted_df


def clean_filename(filepath: str) -> str:
    """
    Clean up a filepath to extract just the filename without extension.

    Args:
        filepath (str): The full file path.

    Returns:
        str: The cleaned filename.
    """
    # Get the base filename without the directory
    base = os.path.basename(filepath)
    
    # Remove the file extension
    name_without_extension = os.path.splitext(base)[0]
    
    return name_without_extension

def main(directory: str):
    """
    Main function to process all MOL2 files in a directory.

    Args:
        directory (str): Path to the directory containing MOL2 files.
    """
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory", file=sys.stderr)
        sys.exit(1)

    df = process_directory(directory)

    df["File"] = [clean_filename(f) for f in df['File']]
    
    # Print summary
    print(f"\nProcessed {df['File'].nunique()} MOL2 files")
    print(f"Found {len(df)} carboxylic acid groups in total")

    # Display the first few rows of the DataFrame
    print("\nFirst few rows of the results:")
    print(df.head())
    df.to_csv("carboxylic_acids_with_charges.csv", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_directory>", file=sys.stderr)
        sys.exit(1)
    
    main(sys.argv[1])