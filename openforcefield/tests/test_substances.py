from functools import partial
import os

from openforcefield import Mixture

def test_neat_liquid():
    liquid = Mixture()
    liquid.addComponent('water')

def test_binary_mixture():
    binary_mixture = Mixture()
    binary_mixture.addComponent('water', mole_fraction=0.2)
    binary_mixture.addComponent('methanol') # assumed to be rest of mixture if no mole_fraction specified

def test_ternary_mixture():
    ternary_mixture = Mixture()
    binary_mixture.addComponent('ethanol', mole_fraction=0.2)
    binary_mixture.addComponent('methanol', mole_fraction=0.2)
    ternary_mixture.addComponent('water')

def test_infinite_dilution():
    infinite_dilution = Mixture()
    infinite_dilution.addComponent('phenol', mole_fraction=0.0) # infinite dilution
    infinite_dilution.addComponent('water')
