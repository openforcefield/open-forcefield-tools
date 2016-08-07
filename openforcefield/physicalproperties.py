#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Physical properties.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>

TODO
----

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

    """
    pass

#=============================================================================================
# THERMOML PHYSICAL PROPERTIES
#=============================================================================================

class MassDensity(PhysicalProperty):
    """Mass density

    """
    pass

class ExcessMolarEnthalpy(PhysicalProperty):
    """Excess molar ExcessMolarEnthalpy

    """
    pass

class HeatCapacity(PhysicalProperty):
    """

    """

#=============================================================================================
# SIMULATION PHYSICAL PROPERTIES OF ISOLATED MOLECULES (USED FOR TESTING ONLY)
#=============================================================================================

class MeanPotentialEnergy(PhysicalProperty):
    """Mean potential energy of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    pass

class MeanPotentialEnergy(PhysicalProperty):
    """Mean potential energy of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    pass

class BondMoment(PhysicalProperty):
    """Specified moment of a specified bond length of an isolated molecule

    The `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified bond of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    pass

class AngleMoment(PhysicalProperty):
    """Specified moment of a specified angle of an isolated molecule

    The `n`th moment (mean `E[x]` if n==1, or central moment `E[(x-E[x])^n]` if n >= 2) of a specified angle of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    pass

class TorsionMoment(PhysicalProperty):
    """Specified circular moment of a specified torsion angle of an isolated molecule

    The `n`th circular moment (`E[e^{i*n*phi}]`) of a specified torsion of an isolated molecule

    Warning
    -------
    This is a simulation-derived property of an isolated molecule, and not a real measurable physical property.

    """
    pass
