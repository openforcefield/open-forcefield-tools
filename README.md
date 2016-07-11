# Open Forcefield Tools

Tools for open forcefield development.

## API documentation

### Physical property measurements

[**Physical property measurements**](https://en.wikipedia.org/wiki/Physical_property) are measured properties of a substance that provide some information about the physical parameters that define the interactions within the substance.

A physical property is defined by a combination of:
* The thermodynamic state the substance that the measurement was performed on, including temperature, pressure, composition, and phase.
* The type of physical property that was was measured (mass density, heat capacity, etc.)
* The experimental approach that was used to obtain the measurement.

An example:

* The thermodynamic state is 298 kelvin, 1 atmosphere, with a mixture of 0.8 mole fraction ethanol and and 0.2 mole fraction water in liquid phase  
* The physical property is the mass density.
* The method that used performed was the vibrating tube method.

In the API, different physical properties such as `MassDensity` are subclasses of the `PhysicalProperty` class.

A `PhysicalProperty` class has properties:
* `.ThermodynamicState`: A `ThermodynamicState` object
* `.MeasurementMethod`: A `MeasurementMethod` object 
* `.value`: the value, a float with units
* `.uncertainty`: the statistical uncertainty, a float with units 
* `.reference` - the literature reference (if present) for the measurement
* `.DOI` - the literature reference DOI (if available) for the measurement
*. `model` - the model that is used to generate the data (if computation)

The value, uncertainty, reference, and DOI do not necessarily need to be defined for a dataset in order for property calculations to be performed.
 
If this physical property is a computation, we could additionally have a 'Model' object, which would include the force field parameters and functional form, but do not need to be defined.

#### Thermodynamic states

A `ThermodynamicState` specifies a combination of thermodynamic parameters (e.g. temperature, pressure) at which a measurement is performed, the phase, as well as the composition.

The `ThermodynamicState` class has properties:
* `.T`: the temperature of the system, a float with simtk units
* `.P`: the pressure of the system, a float with simtk units
* `.composition`: a `Composition` class describing the chemical composition of the system
* `.phase`: the phase, for now a string but will likely be an object eventually.

For interfaces, we many need to add members such as surface tension, etc.

```python
from simtk import unit
thermodynamic_state = ThermodynamicState(pressure=500*unit.kilopascals, temperature=298.15*unit.kelvin, composition=Composition, phase=`liquid`)
```

We use the `simtk.unit` unit system from [OpenMM](http://openmm.org) for units. 

For now, `phase` is a string, and can be `None`, as it does not necessarily need to be used if
carrying out a simulation where there is only one stable equilibrium
state. `phase` can be 'IsolatedMolecule' for simulation properties of
individual molecules.  `liquid` and `gas` are other possible phases.  ThermoML has additional phase descriptions which may be necessary as different systems are investigated (for example, unit cell size and symmetry group for solids).  At some point, it may need to be turned into an object, ideally paralleling [ThermoML](http://trc.nist.gov/ThermoMLRecommendations.pdf) definitions.

The `Composition` class specifies the composition of the system, such as [0.4 mole fraction ethanol, 0.6 mole fraction water].  

We use the concept of a substance's `Composition` throughout, where we
use a class rather than a specific list of molecule fractions or
composition to have more flexibility in how it is specified.

The basic method of the `Composition` class is `addComponent()`, which takes an chemical name as an argument, and an optional mole fraction argument.

A simple liquid has only one component in the `Composition`.
```python
liquid = Composition()
liquid.addComponent('water')
```

A binary mixture has two components:
```python
binary_mixture = Composition()
binary_mixture.addComponent('water', mole_fraction=0.2)
binary_mixture.addComponent('methanol') # assumed to be rest of mixture if no mole_fraction specified
```

A ternary mixture has three components:
```python
ternary_mixture = Composition()
binary_mixture.addComponent('ethanol', mole_fraction=0.2)
binary_mixture.addComponent('methanol', mole_fraction=0.2)
ternary_mixture.addComponent('water')
```
`Composition` is specified in terms of mole fractions of each
  component. If a simulation is performed, some decision will need to
  be performed to decide how large the simulation should be in order
  to get a correct value. Possibly the simulation could also perform
  an extrapolation to larger systems if properly validated.

The infinite dilution of one solute within a solvent or mixture is also specified by setting the `Composition`, where the solute has zero mole fraction:

```python
infinite_dilution = Composition()
infinite_dilution.addComponent('phenol', mole_fraction=0.0) # infinite dilution
infinite_dilution.addComponent('water')
```

i.e., if a component exists, it must have at least one molcule.  Another alternative would be

```python
infinite_dilution = Composition()
infinite_dilution.addComponent('phenol', mole_fraction = 'infinite dilution') # one molecule
infinite_dilution.addComponent('water')  # N-1 water molecules, where N is determined by the simulation.
```

Or something like that.

* Previously, we used the concept of `Mixture` which was an class describing the system and how it was made.  Now, we use the concept of `Composition` which describes how it is made. 

#### Measurement methods

A `MeasurementMethod` subclass has information specific to the particular method used to measure a property (such as experimental uncertainty guidance).  The `MeasurementMethod` object can potentially closely parallel the [ThermoML](http://trc.nist.gov/ThermoMLRecommendations.pdf) specification, which has significant infrastructure for defining the measurement approach.

Some examples:
* `FlowCalorimetry` for `HeatCapacity` or `ExcessMolarEnthalpy`
* `VibratingTubeMethod` for `MassDensity`
* `MolecularSimulation`, where `MolecularSimulation` is an class which contains all the necessary information to generate the physical property given a model (exact enumeration of information TBD). It does not include the parameters for the model.

#### Physical property measurements

* Note: `MeasuredPhysicalProperty` and `ComputedPhysicalProperty` are being merged together into `PhysicalProperty`

An example of creating a `PhysicalProperty`:

```python
# Define mixture
mixture = Composition()
mixture.addComponent('water', mole_fraction=0.2)
mixture.addComponent('methanol')
# Define thermodynamic state
thermodynamic_state = ThermodynamicState(pressure=500*unit.kilopascals, temperature=298.15*unit.kelvin, composition=mixture)
# Define measurement
measurement = ExcessMolarEnthalpy(thermodynamic_state, value=83.3863244*unit.kilojoules_per_mole, uncertainty=0.1220794866*unit.kilojoules_per_mole)
```
The various properties are all subclasses of `PhysicalProperty` and generally follow the `<ePropName/>` ThermoML tag names.
Some examples:
* `MassDensity` - mass density
* `ExcessMolarEnthalpy` - excess partial apparent molar enthalpy
* `HeatCapacity` - molar heat capacity at constant pressure

A [roadmap of physical properties to be implemented](https://github.com/open-forcefield-group/open-forcefield-tools/wiki/Physical-Properties-for-Calculation) is available.
Please raise an issue if your physical property of interest is not listed!

### Physical property datasets

A `PhysicalPropertyDataset` is a collection of `PhysicalProperty` objects that are related in some way.
```python
dataset = PhysicalPropertyDataset([measurement1, measurement2])
```
The dataset is iterable:
```python
dataset = PhysicalPropertyDataset([measurement1, measurement2])
for measurement in dataset:
    print measurement.value
```
and has accessors to retrieve DOIs and references associated with measurements in the dataset:
```python
# Print the DOIs associated with this dataset
print(dataset.DOIs)
# Print the references associated with this dataset
print(dataset.references)
```

### ThermoML datasets

Here, we use `ThermoMLDataset` objects to access datasets stored in the IUPAC-standard [ThermoML](http://trc.nist.gov/ThermoMLRecommendations.pdf) format, a format for specifying thermodynamic properties in XML format. Essentially, `ThermoMLDataset` provides access to ThermoML format datasets.
 
Many datasets from ThermoML itself are of interest, so direct access to the [NIST ThermoML Archive](http://trc.nist.gov/ThermoML.html) is supported for obtaining physical property measurements in this format.

For example, to retrieve [this ThermoML dataset](http://trc.boulder.nist.gov/ThermoML/10.1016/j.jct.2005.03.012) that accompanies [this paper](http://www.sciencedirect.com/science/article/pii/S0021961405000741), we can simply use the DOI `10.1016/j.jct.2005.03.012` as a key for creating a `PhysicalPropertyDataset` subclassed object from the ThermoML Archive:
```python
dataset = ThermoMLDataset('10.1016/j.jct.2005.03.012')
```
You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:
```python
thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
dataset = ThermoMLDataset(thermoml_keys)
```
You can see which DOIs contribute to the current `ThermoMLDataset` with the convenience functions:
```python
thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
dataset = ThermoMLDataset(thermoml_keys)
```
### Estimating properties

The `PropertyEstimator` class creates objects that handle property estimation of all of the properties specified in a dataset, given a set or sets of descriptions of the system.
The implementation will isolate the user from whatever backend (local machine, HPC cluster, XSEDE resources, Amazon EC2) is being used to compute the properties, as well as whether new simulations are being launched and analyzed or existing simulation data is being reweighted.
Different backends will take different optional arguments, but here is an example that will launch and use 10 worker processes on a cluster:
```python
estimator = PropertyEstimator(nworkers=10) # NOTE: multiple backends will be supported in the future
computed_properties = estimator.computeProperties(dataset, model_set)
```
Here, `dataset` is a set of `PhysicalProperties`, and `model_set` is set of a `SMIRFFParameterSet` used to define the physical model of the systems described in the dataset.  There could be one or many elements in the `model_set`. They could differ in force field functional form, or simply parameters.

Optionally, this can take a `MolecularSimulation` object (that can specify how to do the simulation, but this can also be done by default.

`PropertyEstimator.computeProperties(...)` returns a list of `PhysicalProperty` objects. As defined previously, they have:
* `.value` - the computed property value, with appropriate units
* `.uncertainty` - the statistical uncertainty in the computed property
* `.MeasurementMethod` - an object giving the MeasurementMethod (in this case, simulation approach) used to compute this property
* `.Model` - a reference to the parameters used to generate the model. In this case, the model is included, since the properties were generated with the model. Currently, the models will only differ by parameter set, but later the file formats can be updated to allow them to differ in functional form as well.

This API can be extended in the future to provide access to the simulation data used to estimate the property, such as
```python
# Attach to my compute and storage resources
estimator = PropertyEstimator(...)
# Estimate some properties
computed_properties = estimator.computeProperties(dataset, models)
# Get statistics about simulation data that went into each property
for property in computed_properties:
   # Get statistics about simulation data that was reweighted to estimate this property
   for simulation in property.simulations:
      print('The simulation was %.3f ns long' % (simulation.length / unit.nanoseconds))
      print('The simulation was run at %.1f K and %.1f atm' % (simulation.thermodynamic_state.temperature / unit.kelvin, simulation.thermodynamic_state.pressure / unit.atmospheres))
      # Get the ModelSet that was used for this simulation
      parameters = simulation.model
      # what else do you want...?
```

In future, we will want to use a parallel key/value database like [cassandra](http://cassandra.apache.org) to store simulations, along with a distributed task management system like [celery](http://www.celeryproject.org) with [redis](https://www.google.com/search?client=safari&rls=en&q=redis&ie=UTF-8&oe=UTF-8).

## API Usage Examples

### Using the high-level API

In this example, datasets are retrieved from the ThermoML and filtered to retain certain properties.
The corresponding properties for a given parameter set filename are then computed for a SMIRFF parameter set and printed.
```python
# Define the input datasets from ThermoML
thermoml_keys = ['10.1016/j.jct.2005.03.012', ...]
dataset = ThermoMLDataset(thermoml_keys)
# Filter the dataset to include only molar heat capacities measured between 280-350 K
dataset.filter(ePropName='Excess molar enthalpy (molar enthalpy of mixing), kJ/mol') # filter to retain only this property name
dataset.filter(VariableType='eTemperature', min=280*unit.kelvin, max=350*kelvin) # retain only measurements with `eTemperature` in specified range
# Load an initial model set
model_set = [ SMIRFFParameterSet('smarty-initial.xml') ]
# Compute physical properties for these measurements
estimator = PropertyEstimator(nworkers=10) # NOTE: multiple backends will be supported in the future
computed_properties = estimator.computeProperties(dataset, model_set)
# Write out statistics about errors in computed properties
for (computed, measured) in (computed_properties, dataset):
   property_unit = measured.value.unit
   print('%24s : experiment %8.3f (%.3f) | calculated %8.3f (%.3f) %s' % (measured.value / property_unit, measured.uncertainty / property_unit, computed.value / property_unit, computed.uncertainty / property_unit, str(property_unit))
```

### Using the low-level API

TBD
