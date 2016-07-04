#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Physical property estimation API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>

TODO
----
* Implement methods

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
# COMPUTED PROPERTY RESULT
#=============================================================================================

class PropertyComputationSimulation(object):
    """Container for information about a simulation that was run.

    Properties
    ----------
    length : simtk.unit.Quantity with units compatible with nanoseconds
        The length of the simulation.
    thermodynamic_state : ThermodynamicState
        The thermodynamic state at which the simulation was run.
    substance : Substance
        The substance that was simulated.
    system : simtk.openmm.System
        The system that was simulated.

    """
    def __init__(self):
        pass

class ComputedProperty(object):
    """Computed physical property result.

    Properties
    ----------
    value : simtk.unit.Quantity (possibly wrapping a numpy array)
        The estimated (computed) value of the property
    uncertainty : simtk.unit.Quantity with same dimension and units as 'value'
        The estimated uncertainty (standard error) of the computed 'value'
    parameters : ParameterSet
        The parameter set used to compute the provided property
    simulations : set of Simulation objects
        The simulations that contributed to this property estimate


    """
    def __init__(self):
        self.simulations = set()

#=============================================================================================
# PROPERTY ESTIMATOR
#=============================================================================================

class PropertyEstimator(object):
    """Physical property estimation interface.

    This is a generic interface.
    Multiple backends will be supported in the future, and computation is not guaranteed to happen anywhere specific.
    Intermediate files are not guaranteed to be stored; intermediate simulation data may be generated as needed.

    Examples
    --------
    >>> estimator = PropertyEstimator(nworkers=10) # NOTE: multiple backends will be supported in the future
    >>> computed_properties = estimator.computeProperties(dataset, parameter_sets)

    """
    def __init__(self, **kwargs):
        """Create a physical property estimator interface.

        """
        pass

    def computeProperties(self, dataset, parameter_sets):
        """Compute physical properties for the specified dataset given one or more parameter sets.

        Parameters
        ----------
        dataset : PhysicalPropertyDataset
            The dataset for which physical properties are to be computed.
        parameter_sets : ParameterSet or iterable of Parameterset
            Parameter set(s) for which physical properties are to be computed.

        Returns
        -------
        properties : list of ComputedProperty
            The computed

        """
        computed_properties =
        return computed_properties
