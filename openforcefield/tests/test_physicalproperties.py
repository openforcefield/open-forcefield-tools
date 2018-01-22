"""
Tests for PhysicalProperty objects.

"""

from simtk import unit
from openforcefield.thermodynamics import ThermodynamicState
from openforcefield.substances import IsolatedMolecule
from openforcefield.physicalproperties import MeanPotentialEnergy, BondMoment, AngleMoment, TorsionMoment

def test_simulation_properties():
    smiles = 'CCCC'
    molecule = IsolatedMolecule(smiles=smiles)
    thermodynamic_state = ThermodynamicState(temperature=300*unit.kelvin)
    mean_potential = MeanPotentialEnergy(molecule, thermodynamic_state, value=124.4*unit.kilojoules_per_mole, uncertainty=14.5*unit.kilojoules_per_mole)
    bond_average = BondMoment(molecule, thermodynamic_state, value=1.12*unit.angstroms, uncertainty=0.02*unit.angstroms, moment=1, smirks='[#6:1]-[#6:2]')
    bond_variance = BondMoment(molecule, thermodynamic_state, value=0.05*unit.angstroms**2, uncertainty=0.02*unit.angstroms**2, moment=2, smirks='[#6:1]-[#6:2]')
    angle_average = AngleMoment(molecule, thermodynamic_state, value=20*unit.degrees, uncertainty=0.05*unit.degrees, moment=1, smirks='[#6:1]-[#6:2]-[#6:3]')
    torsion_moment = TorsionMoment(molecule, thermodynamic_state, value=(0.5, 0.3), uncertainty=(0.05, 0.03), moment=1, smirks='[#6:1]-[#6:2]-[#6:3]-[#6:4]')
