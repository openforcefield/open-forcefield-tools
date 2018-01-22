#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Substances API.

Authors
-------
* John D. Chodera <john.chodera@choderalab.org>

TODO
----
* Incorporate Ken Kroenlein's suggestion that infinite dilution should be treated differently than zero mole fraction.
  Perhaps addInfiniteDilution or addImpurity?
* Add methods that construct real System and Topology objects for a specified system size, following the Mobley SolvationToolkit:
  https://github.com/MobleyLab/SolvationToolkit

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
from openeye import oechem, oeiupac, oeomega

from distutils.spawn import find_executable
from collections import OrderedDict

#=============================================================================================
# SUBSTANCE
#=============================================================================================

class Substance(object):
    """A substance.
    """
    pass

#=============================================================================================
# ISOLATED MOLECULE
#=============================================================================================

class IsolatedMolecule(Substance):
    """An isolated molecule.

    This is used only for testing.

    Examples
    --------

    Specify a molecule using SMILES

    >>> ethane = IsolatedMolecule(smiles='CC')

    Specify a molecule using IUPAC

    >>> ethane = IsolatedMolecule(iupac='ethane')

    """
    def __init__(self, smiles=None, iupac=None):
        """Create an isolated molecule.

        Parameters
        ----------
        smiles : str, optional, default=None
            An isomeric SMILES string specifying the isolated molecule.
        iupac : str, optional, default=None
            The IUPAC name of the isolated molecule.

        One of `smiles` or `iupac` must be specified
        """
        if (smiles is None) and (iupac is None):
            raise Exception("Either 'smiles' or 'iupac' must be specified.")
        if (smiles is not None) and (iupac is not None):
            raise Exception("Only one of 'smiles' and 'iupac' can be specified.")
        self.smiles = smiles
        self.iupac = iupac

#=============================================================================================
# MIXTURE
#=============================================================================================

class Mixture(Substance):
    """A liquid or gas mixture.

    Properties
    ----------
    components : dict
        components[iupac_name] is the mole fraction of the specified component

    Examples
    --------

    A neat liquid has only one component:

    >>> liquid = Mixture()
    >>> liquid.addComponent('water')

    A binary mixture has two components:

    >>> binary_mixture = Mixture()
    >>> binary_mixture.addComponent('water', mole_fraction=0.2)
    >>> binary_mixture.addComponent('methanol') # assumed to be rest of mixture if no mole_fraction specified

    A ternary mixture has three components:

    >>> ternary_mixture = Mixture()
    >>> binary_mixture.addComponent('ethanol', mole_fraction=0.2)
    >>> binary_mixture.addComponent('methanol', mole_fraction=0.2)
    >>> ternary_mixture.addComponent('water')

    The infinite dilution of one solute within a solvent or mixture is also specified as a `Mixture`, where the solute has zero mole fraction:

    >>> infinite_dilution = Mixture()
    >>> infinite_dilution.addComponent('phenol', mole_fraction=0.0) # infinite dilution
    >>> infinite_dilution.addComponent('water')

    """

    class Component(object):
        _cached_molecules = dict() # store cached molecules by IUPAC name
        def _createMolecule(self, iupac_name):
            """Create molecule from IUPAC name.

            Best practices for creating initial coordinates can be applied here.

            Parameters
            ----------
            iupac_name : str
                IUPAC name

            Returns
            -------
            molecule : OEMol
                OEMol with 3D coordinates, but no charges

            """
            # Check cache
            if iupac_name in self._cached_molecules:
                return copy.deepcopy(self._cached_molecules[iupac_name])

            # Create molecule from IUPAC name.
            molecule = oechem.OEMol()
            if not oeiupac.OEParseIUPACName(molecule, iupac_name):
                raise ValueError("The supplied IUPAC name '%s' could not be parsed." % iupac_name)

            # Set molecule name
            molecule.SetTitle(iupac_name)

            # Normalize molecule
            oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)
            oechem.OEAddExplicitHydrogens(molecule)
            oechem.OETriposAtomNames(molecule)

            # Create configuration
            omega = oeomega.OEOmega()
            omega.SetMaxConfs(1)
            omega.SetIncludeInput(False)
            omega.SetCanonOrder(False)
            omega.SetSampleHydrogens(True)
            omega.SetStrictStereo(True)
            omega.SetStrictAtomTypes(False)
            status = omega(molecule)
            if not status:
                raise(RuntimeError("omega returned error code %d" % status))

            # Cache molecule
            self._cached_molecules[iupac_name] = molecule

            return molecule

        def __init__(self, iupac_name, mole_fraction=0.0, impurity=False):
            """Create a mixture component.

            Parameters
            ----------
            iupac_name : str
                IUPAC name of comonent
            mole_fraction : float, optional, default=0.0
                Mole fraction (must be zero if impurity=True)
            impurity : bool, optional, default=False
                If True, this is an impurity

            """
            # TODO: Check that IUPAC name is not stereochemically ambiguous.

            if not ((mole_fraction >= 0.0) and (mole_fraction <= 1.0)):
                raise Exception('Mole fraction must be positive; specified %f.' % mole_fraction)
            if impurity and (mole_fraction != 0.0):
                raise Exception('Mole fraction must be 0.0 for impurities. Specified mole fraction of %f' % mole_fraction)

            self.iupac_name = iupac_name
            self.mole_fraction = mole_fraction
            self.impurity = impurity
            self.molecule = self._createMolecule(iupac_name)

    def __init__(self):
        """Create a Mixture.
        """
        self.components = list()

    @property
    def total_mole_fraction(self):
        """Compute the total mole fraction.
        """
        return sum([ component.mole_fraction for component in self.components ])

    @property
    def n_components(self):
        return len(self.components)

    @property
    def n_impurities(self):
        return sum([ 1 for component in self.components if component.impurity==True ])

    def addComponent(self, iupac_name, mole_fraction=None, impurity=False):
        """Add a component to the mixture.

        Parameters
        ----------
        iupac_name : str
            IUPAC name of the component.
        mole_fraction : float, optional, default=None
            If specified, the mole fraction of this component.
            If not specified, this will be the last or only component of the mixture.
        impurity : bool, optional, default=False
            If True, the component represents an impurity (with zero mole_fraction).

        """
        if impurity and (mole_fraction is None):
            mole_fraction = 0.0

        if mole_fraction is None:
            mole_fraction = 1.0 - self.total_mole_fraction

        if (self.total_mole_fraction + mole_fraction) > 1.0:
            raise Exception("Total mole fraction would exceed unity (%f); specified %f." % (total_mole_fraction, mole_fraction))

        component = Mixture.Component(iupac_name=iupac_name, mole_fraction=mole_fraction, impurity=impurity)
        self.components.append(component)

    def getComponent(self, iupac_name):
        """Retrieve component by IUPAC name.

        Parameters
        ----------
        iupac_name : str
            The name of the component to retrieve

        """
        for component in self.components:
            if component.iupac_name == iupac_name:
                return component
        raise Exception("No component with IUPAC name '%s' found." % iupac_name)

    def build(self, nmolecules=1000, mass_density=None):
        """Build an instance of this mixture.

        Parameters
        ----------
        nmolecules : int, optional, default=True
            The number of molecules in the system to be created.
        mass_density : float, optional, default=True
            If provided, will aid in the selecting an initial box size.

        Returns
        -------
        topology : simtk.openmm.Topology
            Topology object describing the system.
        molecules : list of oechem.OEMol
            The molecules describing the system (not guaranteed to have charges or 3D coordinates).
            These are copies of molecules.
        positions : simtk.unit.Quantity wrapping [natoms,3] numpy array with units compatible with angstroms
            Positions of all atoms in the system.

        Notes
        -----
        The number of molecules of each component need not be deterministic.
        Repeated calls may generate different numbers of molecules, orders, and positions.
        Impurities will have exactly one molecule per impurity.

        """

        # Create deep copies of molecules.
        molecules = [ copy.deepcopy(component.molecule) for component in self.components ]

        # Determine how many molecules of each type will be present in the system.
        mole_fractions = np.array([ component.mole_fraction for component in self.components ])
        n_copies = np.random.multinomial(nmolecules - self.n_impurities, pvals=mole_fractions)

        # Each impurity must have exactly one molecule
        for (index, component) in enumerate(self.components):
            if component.impurity:
                n_copies[index] = 1

        # Create packed box
        from openforcefield.packmol import pack_box
        [topology, positions] = pack_box(molecules, n_copies, mass_density=mass_density)

        return [topology, molecules, positions]
