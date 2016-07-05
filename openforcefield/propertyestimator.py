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

class ComputedPropertySet(list):
    """Set of computed physical properties.

    """
    def __init__(self):
        pass

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

    def computeProperties(self, dataset, parametersets, target_relative_uncertainty=0.1):
        """Compute physical properties for the specified dataset given one or more parameter sets.

        Parameters
        ----------
        dataset : PhysicalPropertyDataset
            The dataset for which physical properties are to be computed.
        parametersets : ParameterSet or iterable of Parameterset
            Parameter set(s) for which physical properties are to be computed.
        target_uncertainty : float, optional, default=0.1
            Target computational uncertainty in the computed property, relative to experimental uncertainty.

        Returns
        -------
        properties : ComputedPropertySet object or list of ComputedPropertySet objects
            The computed physical properties.

        """
        # Attempt to estimate all simulated properties by reweighting
        simulations_to_run = list() # list of simulations that need to be rerun
        computed_properties_sets = list()
        for parameters in parameter_sets:
            computed_property_set = ComputedPropertySet()
            for measured_property in dataset:
                # Estimate property via reweighting.
                computed_property = self._estimateProperty(measured_property)
                if (computed_property is None) or (computed_property.uncertainty > target_relative_uncertainty * measured_property.uncertainty):
                    # Uncertainty threshold exceeded; queue for simulation
                    simulation = Simulation(thermodynamic_state=measured_property.thermodynamic_state, composition=thermodynamic_state.composition)
                    simulations_to_run.append( (simulation, measured_property, target_relative_uncertainty) )

                computed_property = ComputedProperty()
                computed_property_set.append(computed_property)
            computed_properties_sets.append(computed_property_set)

        # Run queued simulations
        # TODO: Parallelize
        for (simulation, target_property, target_uncertainty) in simulations_to_run:
            simulation.run(target_property=target_property, target_uncertainty=target_uncertainty)

        # Return computed physical property dataset(s)
        return computed_properties_sets
