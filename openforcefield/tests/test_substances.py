from functools import partial
import os

from openforcefield.substances import IsolatedMolecule, Mixture

def test_isolated_molecule():
    ethane = IsolatedMolecule(smiles='CC')
    ethane = IsolatedMolecule(iupac='ethane')

def test_component():
    water = Mixture.Component('water', mole_fraction=1.0)
    methanol = Mixture.Component('methanol', mole_fraction=0.2)

def test_neat_liquid():
    mixture = Mixture()
    mixture.addComponent('water')
    assert (mixture.getComponent('water').mole_fraction == 1.0)
    assert (mixture.n_components == 1)
    [topology, molecules, positions] = mixture.build(nmolecules=1000)

def test_binary_mixture():
    mixture = Mixture()
    mixture.addComponent('water', mole_fraction=0.2)
    mixture.addComponent('methanol')
    assert (mixture.getComponent('water').mole_fraction == 0.2)
    assert (mixture.getComponent('methanol').mole_fraction == 0.8)
    assert (mixture.n_components == 2)
    [topology, molecules, positions] = mixture.build(nmolecules=1000)

def test_ternary_mixture():
    mixture = Mixture()
    mixture.addComponent('ethanol', mole_fraction=0.2)
    mixture.addComponent('methanol', mole_fraction=0.2)
    mixture.addComponent('water')
    assert (mixture.getComponent('ethanol').mole_fraction == 0.2)
    assert (mixture.getComponent('methanol').mole_fraction == 0.2)
    assert (mixture.getComponent('water').mole_fraction == 0.6)
    assert (mixture.n_components == 3)
    [topology, molecules, positions] = mixture.build(nmolecules=1000)

def test_infinite_dilution():
    mixture = Mixture()
    mixture.addComponent('phenol', impurity=True)
    mixture.addComponent('water')
    assert (mixture.getComponent('phenol').mole_fraction == 0.0)
    assert (mixture.getComponent('phenol').impurity == True)
    assert (mixture.getComponent('water').impurity == False)
    assert (mixture.n_components == 2)
    assert (mixture.n_impurities == 1)
    [topology, molecules, positions] = mixture.build(nmolecules=1000)

def test_component_iteration():
    mixture = Mixture()
    mixture.addComponent('ethanol', mole_fraction=0.2)
    mixture.addComponent('methanol', mole_fraction=0.2)
    mixture.addComponent('water')

    for component in mixture.components:
        iupac_name = component.iupac_name
        mole_fraction = component.mole_fraction
        impurity = component.impurity

def test_create_instance():
    mixture = Mixture()
    mixture.addComponent('ethanol', mole_fraction=0.2)
    mixture.addComponent('methanol', mole_fraction=0.2)
    mixture.addComponent('water')

    for nmolecules in [10, 100, 1000, 10000]:
        [topology, molecules, positions] = mixture.build(nmolecules=1000)
        f = partial(mixture.build, nmolecules)
        f.description = "Testing Mixture.createInstance for %d molecules" % nmolecules
        yield f
