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
    """A liquid mixture.

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
    def __init__(self):
        """Create a Mixture.
        """
        self.components = dict()

    def _totalMoleFraction(self):
        """Compute the total mole fraction.
        """
        total_mole_fraction = 0.0
        for component in self.components:
            total_mole_fraction += self.components[component]
        return total_mole_fraction

    def addComponent(self, component, mole_fraction=None):
        """Add a component to the mixture.

        Parameters
        ----------
        component : str
            IUPAC name of the component.
        mole_fraction : float, optional, default=None
            If specified, the mole fraction of this component.
            If not specified, this will be the last or onlycomponent of the mixture.

        """
        total_mole_fraction = self._total_mole_fraction
        if mole_fraction is not None:
            if (total_mole_fraction + mole_fraction) > 1.0:
                raise Exception("Total mole fraction would exceed unity (%f)" % total_mole_fraction)
            self.components[component] = mole_fraction
        else:
            # Utilize rest of mole fraction
            self.components[component] = 1.0 - total_mole_fraction
