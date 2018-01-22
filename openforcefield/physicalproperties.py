#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Physical properties.

Right now, these classes serve mainly to identify the type of measurement carried out.
Eventually, they will also hold utility functions for computing physical properties from simulations.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>

TODO
----
* Consider general MapReduce framework for each property's analysis scheme.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import sys
import time
import copy
import numpy as np

from simtk import openmm, unit

#=============================================================================================
# PHYSICAL PROPERTY BASE CLASS
#=============================================================================================

class PhysicalProperty(object):
    """A physical property

    Properties
    ----------


    """
    pass

#=============================================================================================
# ESTIMATION OF EXPERIMENTAL OBSERVABLES
#=============================================================================================

class ObservableEstimationStrategy(object):
    """
    """
    def __init__(self, substance, topology, system):
        """
        """
        pass

    def map(self, substance, topology, system, state):
        """
        """
        pass

    def reduce(self, map_result, thermodynamic_state):
        """
        """
        pass

#=============================================================================================
# THERMOML PHYSICAL PROPERTIES
#=============================================================================================

class MassDensity(PhysicalProperty):
    """Mass density.

    """

    def _computeMass(self, topology):
        """Compute total mass of the system.

        TODO
        ----
        * Use memoization to speed this up?
        """
        mass = 0.0 * unit.amu
        for atom in topology.atoms():
            mass += atom.element.mass

        return mass

    def computeObservable(self, substance, topology, system, state):
        """Compute the relevant observable

        """
        mass = self._computeMass(topology)
        volume = state.getPeriodicBoxVolume()
        density = mass/volume
        return density

    pass

class ExcessMolarEnthalpy(PhysicalProperty):
    """Excess molar enthalpy.

    """
    pass

class HeatCapacity(PhysicalProperty):
    """Heat capacity.

    """
    pass

class StaticDielectricConstant(PhysicalProperty):
    """Static dielectric constant.

    """
    pass

#=============================================================================================
# SIMULATION PHYSICAL PROPERTIES OF ISOLATED MOLECULES (USED FOR TESTING ONLY)
#=============================================================================================

class MeanPotentialEnergy(PhysicalProperty):
    """Mean potential energy of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    def __init__(self, molecule, thermodynamic_state, value, uncertainty):
        self.thermodynamic_state = thermodynamic_state
        self.value = value
        self.uncertainty = uncertainty

class BondMoment(PhysicalProperty):
    """Specified moment of a specified bond length of an isolated molecule

    The `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified bond of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    def __init__(self, molecule, thermodynamic_state, value, uncertainty, moment, smirks):
        self.thermodynamic_state = thermodynamic_state
        self.value = value
        self.uncertainty = uncertainty

class AngleMoment(PhysicalProperty):
    """Specified moment of a specified angle of an isolated molecule

    The `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified angle of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    def __init__(self, molecule, thermodynamic_state, value, uncertainty, moment, smirks):
        self.thermodynamic_state = thermodynamic_state
        self.value = value
        self.uncertainty = uncertainty

class TorsionMoment(PhysicalProperty):
    """Specified circular moment of a specified torsion angle of an isolated molecule

    The `n`th circular moment (`E[e^{i*n*phi}]`) of a specified torsion of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    def __init__(self, molecule, thermodynamic_state, value, uncertainty, moment, smirks):
        self.thermodynamic_state = thermodynamic_state
        self.value = value
        self.uncertainty = uncertainty
