#!/bin/env python

"""Utility functions for use with other modules here.

Author: D. Mobley
"""

from smarty import ForceField
from smarty import *
import openeye
from openeye import oechem
from openeye import oeomega
from openeye import oequacpac
import smarty
from simtk import openmm
from simtk.openmm import app
from simtk import unit
import numpy as np
from smarty.utils import get_data_filename
import random
import matplotlib
import pylab as pl
import os


# Define utility function we'll use to get energy of an OpenMM system
def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy (kcal/mol)
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    return energy

def get_minimum_energy(system, topology, positions):
    """
    Return the potential energy after minimization, in kcal/mol.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    topology : OpenMM topology
        Topology for system
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy (kcal/mol)
    """

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    simulation=app.Simulation( topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    return energy

def reformat_oemol_coordinates( oemol ):
    """
    Take an oemol with multiple conformers and return a conformer_positions numpy array which is an Nx3xM array
    where N is the number of atoms, each has three coordinates, and M is the number of conformers.

    Parameters
    ----------
    oemol : OpenEye oemol, multiconformer
        Multi-conformer molecule to work with

    Returns
    ----------
    conformer_positions : numpy array
        Nx3xM array of simtk.unit.Quantity of dimension (natoms,3, M) where M is the number of quantity, with units of length
        """

    #Get positions for use below
    numconfs = oemol.NumConfs()
    coordinates = oemol.GetCoords()
    natoms=len(coordinates)

    conformer_positions = np.zeros([natoms,3, numconfs], np.float32)
    for (idx,conf) in enumerate(oemol.GetConfs()):
        for index in range(natoms):
            coordinates=conf.GetCoords()
            (x,y,z) = coordinates[index]
            conformer_positions[index,0, idx] = x
            conformer_positions[index,1, idx] = y
            conformer_positions[index,2, idx] = z
    conformer_positions = unit.Quantity(conformer_positions, unit.angstroms)
    return conformer_positions

def log_likelihood( z, mu, s_E):
    """Calculate a log likelihood and a likelihood given a calculated value, a set of `measured`/observed values
    (data we are fitting to) and a standard error in those measured values.


    Parameters
    ----------
    z : simtk.unit.Quantity of dimensions N
        Calculated values
    mu : simtk.unit.Quantity of dimensions N
        Mean (i.e. observed/measured values)
    s_E : simtk.unit.Quantity of dimensions N
        Standard error in the mean (uncertainty in observed/measured values)

    Returns
    -------
    loglike : float
        Natural logarithm of likelihood
    like : float
        Likelihood
    """
    # Adapted from https://github.com/shirtsgroup/lj_bayesian/blob/041b896d37f91f4b42cccb2df73af84a9cf5b917/generate_posterior.py#L117


    # Standardize units of input and strip units, if they have units
    try:
        unit_choice = z.unit
        clear_s_E = s_E.value_in_unit(unit_choice)
        clear_mu = mu.value_in_unit(unit_choice)
        clear_z = z.value_in_unit(unit_choice)
    except AttributeError:
        clear_s_E = s_E
        clear_mu = mu
        clear_z = z
    # If we've mixed things with units with things without units, the below will raise another attribute error

    # Compute log likelihood and likelihood
    # Gaussian is (2*pi*s_E)^(-1/2) e^(-(z-mu)^2 /2s_E^2)
    # log Gaussian is -1/2 * (log(2*pi) + log(s_E)) -|z-mu|^2/(2s_E^2)
    # Here we will ignore the constant term because it drops out when taking ratios of probabilities
    loglike = 0.
    like = 0.
    for i, m in enumerate(clear_mu):
        this_s_E = clear_s_E[i]
        this_z = clear_z[i]
        gauss_arg = - ((m-this_z)**2)/(2.*this_s_E**2)
        loglike += gauss_arg - 0.5*np.log(2*np.pi*this_s_E)

    like = np.exp(loglike)

    return loglike, like


