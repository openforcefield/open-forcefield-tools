#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Thermodynamics API.

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
# THERMODYNAMIC STATE
#=============================================================================================

class ThermodynamicState(object):
    """
    Data specifying a thermodynamic state obeying Boltzmann statistics.

    Properties
    ----------
    temperature : simtk.unit.Quantity with units compatible with kelvin
        The external temperature
    pressure : simtk.unit.Quantity with units compatible with atmospheres
        The external pressure

    Examples
    --------
    Specify an NPT state at 298 K and 1 atm pressure.

    >>> state = ThermodynamicState(temperature=298.0*unit.kelvin, pressure=1.0*unit.atmospheres)

    Note that the pressure is only relevant for periodic systems.

    """

    def __init__(self, temperature=None, pressure=None):
        """
        Initialize the thermodynamic state.

        Parameters
        ----------
        temperature : simtk.unit.Quantity compatible with 'kelvin', optional, default=None
           The temperature for a system with constant temperature
        pressure : simtk.unit.Quantity compatible with 'atmospheres', optional, default=None
           The pressure for constant-pressure systems (default: None)

        """
        self.temperature = temperature
        self.pressure = pressure

        # Check that units are compatible.
        # These throw a TypeException if this conversion cannot be done.
        temperature.in_units_of(unit.kelvin)
        pressure.in_units_of(unit.atmospheres)

    def is_compatible_with(self, state):
        """
        Determine whether another state is in the same thermodynamic ensemble (e.g. NVT, NPT).

        Parameters
        ----------
        state : ThermodynamicState
            Thermodynamic state whose compatibility is to be determined

        Returns
        -------
        is_compatible : bool
            True if 'state' is of the same ensemble (e.g. both NVT, both NPT), False otherwise

        Examples
        --------

        Create NVT and NPT states.

        >>> from simtk import unit
        >>> from openmmtools import testsystems
        >>> testsystem = testsystems.LennardJonesCluster()
        >>> [system, positions] = [testsystem.system, testsystem.positions]
        >>> nvt_state = ThermodynamicState(system=system, temperature=100.0*unit.kelvin)
        >>> npt_state = ThermodynamicState(system=system, temperature=100.0*unit.kelvin, pressure=1.0*unit.atmospheres)

        Test compatibility.

        >>> test1 = nvt_state.is_compatible_with(nvt_state)
        >>> test2 = nvt_state.is_compatible_with(npt_state)
        >>> test3 = npt_state.is_compatible_with(nvt_state)
        >>> test4 = npt_state.is_compatible_with(npt_state)

        """

        is_compatible = True

        # Make sure systems have the same number of atoms.
        if ((self.system != None) and (state.system != None)):
            if (self.system.getNumParticles() != state.system.getNumParticles()):
                is_compatible = False

        # Make sure other terms are defined for both states.
        # TODO: Use introspection to get list of parameters?
        for parameter in ['temperature', 'pressure']:
            if (parameter in dir(self)) is not (parameter in dir(state)):
                # parameter is not shared by both states
                is_compatible = False

        return is_compatible

    def __repr__(self):
        """
        Returns a string representation of a state.

        Examples
        --------
        Create an NVT state.

        >>> from simtk import unit
        >>> from openmmtools import testsystems
        >>> testsystem = testsystems.LennardJonesCluster()
        >>> [system, positions] = [testsystem.system, testsystem.positions]
        >>> state = ThermodynamicState(system=system, temperature=100.0*unit.kelvin)

        Return a representation of the state.

        >>> state_string = repr(state)

        """

        r = "ThermodynamicState("
        narguments = 0
        if self.temperature is not None:
            if narguments > 0: r += ", "
            r += "temperature=%s" % repr(self.temperature)
            narguments += 1
        if self.pressure is not None:
            if narguments > 0: r += ", "
            r += ", pressure = %s" % repr(self.pressure)
            narguments += 1
        r += ")"

        return r

    def __str__(self):

        r = "<ThermodynamicState object"
        if self.temperature is not None:
            r += ", temperature = %s" % str(self.temperature)
        if self.pressure is not None:
            r += ", pressure = %s" % str(self.pressure)
        r += ">"

        return r
