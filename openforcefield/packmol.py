"""
Packmol utility routines.

"""
import numpy as np
import shutil
import os
import tempfile
import copy
from distutils.spawn import find_executable

from simtk import unit
from simtk.openmm import app
from openeye import oechem

PACKMOL_PATH = find_executable("packmol")

HEADER_TEMPLATE = """
# Mixture
tolerance %f
filetype pdb
output %s
add_amber_ter
"""

BOX_TEMPLATE = """
structure %s
  number %d
  inside box 0. 0. 0. %f %f %f
end structure
"""

def pack_box(molecules, n_copies, tolerance=2.0, box_size=None, mass_density=None, verbose=False):
    """Run packmol to generate a box containing a mixture of molecules.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the molecules
    tolerance : float, optional, default=2.0
        The mininum spacing between molecules during packing.  In ANGSTROMS!
    box_size : simtk.unit.Quantity in units compatible with angstroms
        The size of the box to generate.
        Default generates boxes that are very large for increased stability.
        May require extra time for energy minimization and equilibration.
    mass_density : simtk.unit.Quantity with units compatible with grams/milliliters, optional, default = 1.0*grams/milliliters
        Target mass density for final system, if available.
    verbose : bool, optional, default=False
        If True, verbose output is written.

    Returns
    -------
    topology : simtk.openmm.Topology
        Topology of the resulting system
    positions : simtk.unit.Quantity wrapped [natoms,3] numpy array with units compatible with angstroms
        Single frame trajectory with mixture box.

    """
    assert len(molecules) == len(n_copies), "Length of 'molecules' and 'n_copies' must be identical"

    # Create PDB files for all components
    pdb_filenames = list()
    for molecule in molecules:
        tmp_filename = tempfile.mktemp(suffix=".pdb")
        pdb_filenames.append(tmp_filename)
        ofs = oechem.oemolostream(tmp_filename)
        oechem.OEWriteMolecule(ofs, copy.deepcopy(molecule))
        ofs.close()

    # Run packmol
    PACKMOL_PATH = find_executable("packmol")
    if PACKMOL_PATH is None:
        raise(IOError("Packmol not found, cannot run pack_box()"))

    output_filename = tempfile.mktemp(suffix=".pdb")

    # Approximate volume to initialize box
    if (box_size is None):
        if (mass_density is not None):
            # Estimate box_size from mass density.
            box_size = approximate_volume_by_density(molecules, n_copies, mass_density=mass_density)
        else:
            # Use vdW radii to estimate box_size
            box_size = approximate_volume(molecules, n_copies)

    header = HEADER_TEMPLATE % (tolerance, output_filename)
    for (pdb_filename, molecule, count) in zip(pdb_filenames, molecules, n_copies):
        header = header + BOX_TEMPLATE % (pdb_filename, count, box_size / unit.angstroms, box_size / unit.angstroms, box_size / unit.angstroms)

    pwd = os.getcwd()

    if verbose: print(header)

    # Write packmol input
    packmol_filename = "packmol_input.txt"
    packmol_filename = tempfile.mktemp(suffix=".txt")
    file_handle = open(packmol_filename, 'w')
    file_handle.write(header)
    file_handle.close()

    os.system("%s < %s" % (PACKMOL_PATH, packmol_filename))

    # Read the resulting PDB file.
    pdbfile = app.PDBFile(output_filename)

    # Extract topology and positions
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()

    return [topology, positions]

def approximate_volume(molecules, n_copies, box_scaleup_factor=2.0):
    """Approximate the appropriate box size based on the number and types of atoms present.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the molecules
    box_scaleup_factor : float, optional, default = 2.0
        Factor by which the estimated box size is increased

    Returns
    -------
    box_size : simtk.unit.Quantity with units compatible with angstroms
        The size of the box to generate.

    Notes
    -----
    By default, boxes are very large for increased stability, and therefore may
    require extra time for energy minimization and equilibration.

    """
    volume = 0.0 * unit.angstrom**3
    for (molecule, number) in zip(molecules, n_copies):
        molecule_volume = 0.0 * unit.angstrom**3
        for atom in molecule.GetAtoms():
            molecule_volume += oechem.OEGetBondiVdWRadius(atom.GetAtomicNum()) * unit.angstrom**3
        volume += molecule_volume * number
    box_edge = volume**(1.0/3.0) * box_scaleup_factor

    return box_edge

def approximate_volume_by_density(molecules, n_copies, mass_density=1.0*unit.grams/unit.milliliters, box_scaleup_factor=1.1):
    """Generate an approximate box size based on the number and molecular weight of molecules present, and a target density for the final solvated mixture. If no density is specified, the target density is assumed to be 1 g/ml.

    Parameters
    ----------
    molecules : list of OEMol
        Molecules in the system (with 3D geometries)
    n_copies : list of int (same length as 'molecules')
        Number of copies of the moleculesixture component.
    box_scaleup_factor : float, optional, default = 1.1
        Factor by which the estimated box size is increased
    mass_density : simtk.unit.Quantity with units compatible with grams/milliliters, optional, default = 1.0*grams/milliliters
        Target mass density for final system, if available.

    Returns
    -------
    box_edge : simtk.unit.Quantity with units compatible with angstroms
        The size (edge length) of the box to generate.

    Notes
    -----
    By default, boxes are only modestly large. This approach has not been extensively tested for stability but has been used in th Mobley lab for perhaps ~100 different systems without substantial problems.

    """
    # Load molecules to get molecular weights
    volume = 0.0 * unit.angstrom**3
    for (molecule, number) in zip(molecules, n_copies):
        molecule_mass = oechem.OECalculateMolecularWeight(mol) * unit.grams/units.mole / unit.AVOGADRO_CONSTANT_NA
        molecule_volume = molecule_mass / mass_density
        volume += molecule_volume * number

    box_edge = volume**(1.0/3.0) * box_scaleup_factor

    return box_edge
