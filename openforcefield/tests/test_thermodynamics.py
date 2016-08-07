from functools import partial
import os

from simtk import unit
from openforcefield.thermodynamics import ThermodynamicState

def test_nvt():
    """Creating an NVT ThermodynamicState"""
    nvt = ThermodynamicState(temperature=300*unit.kelvin)

def test_npt():
    """Creating an NPT ThermodynamicState"""
    npt = ThermodynamicState(temperature=300*unit.kelvin, pressure=1.0*unit.atmosphere)

@raises(TypeException)
def test_wrong_units():
    """Creating a ThermodynamicState with incorrect units"""
    nvt = ThermodynamicState(temperature=300*unit.meters)
