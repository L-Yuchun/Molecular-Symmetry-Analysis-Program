import numpy as np
import argparse
from pymatgen.core import Molecule, Structure
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
from ase import io

def read_structure(file_path):
    """
    Read the molecular or crystal structure from various file formats.

    Parameters:
        file_path (str): Path to the structure file.

    Returns:
        structure (Structure or Molecule): The structure object from pymatgen.
    """
    try:
        # Attempt to read using pymatgen for CIF, POSCAR, and VASP files
        if file_path.lower().endswith((".cif", ".poscar", ".vasp")):
            structure = Structure.from_file(file_path)
        else:
            # Use ASE to read other formats and convert to pymatgen objects
            ase_structure = io.read(file_path)
            positions = np.array(ase_structure.get_positions(), dtype=float)
            symbols = ase_structure.get_chemical_symbols()

            # Check for periodic boundary conditions
            if np.any(ase_structure.get_pbc()):
                lattice = ase_structure.get_cell()
                structure = Structure(lattice, symbols, positions)
            else:
                structure = Molecule(symbols, positions)
        return structure
    except Exception as e:
        print(f"Error reading the structure file: {e}")
        return None

def analyze_symmetry(structure):
    """
    Analyze the point group symmetry of the molecule or crystal.

    Parameters:
        structure (Structure or Molecule): The structure object.

    Returns:
        point_group (str): The point group of the molecule or crystal.
    """
    if isinstance(structure, Molecule):
        # Directly analyze the point group for Molecule objects
        analyzer = PointGroupAnalyzer(structure)
        point_group = analyzer.sch_symbol  # Use sch_symbol for point group
    elif isinstance(structure, Structure):
        # For Structure objects, treat it as a molecule if no PBC
        if not structure.is_periodic:
            molecule = Molecule(structure.species, structure.cart_coords)
            analyzer = PointGroupAnalyzer(molecule)
            point_group = analyzer.sch_symbol  # Use sch_symbol for point group
        else:
            print("The structure is periodic; point group analysis is not applicable.")
            point_group = "Periodic Structure"
    else:
        point_group = "Unknown"
    return point_group

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Analyze molecular symmetry from a structure file.")
    parser.add_argument("file_path", type=str, help="Path to the structure file (e.g., .cif, .pdb, .mol)")
    
    # Parse command-line arguments
    args = parser.parse_args()
    
    # Read the structure from the provided file
    structure = read_structure(args.file_path)
    
    if structure:
        # Analyze the point group symmetry
        point_group = analyze_symmetry(structure)
        # Output the symmetry information
        print(f"The point group of the molecule is: {point_group}")
    else:
        print("Failed to read the structure file.")
