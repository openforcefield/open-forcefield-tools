from functools import partial
import os

from simtk import unit

SIGMA_CUTOFF = 6.0 # raise an Exception when computed properties differ from measured properties by more than this multiple of statistical error

def test_propertyestimator():
    """Test PhysicalPropertyEstimator on synthetic data.

    """
    # Get synthetic dataset and parameters
    from openforcefield import testsystems
    dataset = testsystems.TestDataset()
    parameters = testsystems.TestParameterSet()

    # Create a PropertyEstimator
    from openforcefield import PropertyEstimator
    estimator = PropertyEstimator()

    # Compute properties
    computed_properties = estimator.computeProperties(dataset, parameters)

    # Assert that computed properties are within statistical error
    for (measured_property, computed_property) in zip(dataset, computed_properties):
        error_value = (computed_property.value - measured_property.value)
        error_uncertainty = unit.sqrt(computed_property.uncertainty**2 + measured_property.uncertainty**2)
        relative_error = error_value / error_uncertainty
        if (relative_error > SIGMA_CUTOFF):
            msg  = 'Computed property %s differs from measured property by more than SIGMA_CUTOFF (%f):\n' % (computed_property.name, SIGMA_CUTOFF)
            msg += 'Measured: %12.3f +- %12.3f %s' % (measured_property.value / measured_property.unit, measured_property.uncertainty / measured_property.unit, str(measured_property.unit))
            msg += 'Computed: %12.3f +- %12.3f %s' % (computed_property.value / measured_property.unit, computed_property.uncertainty / measured_property.unit, str(measured_property.unit))
            msg += 'ERROR   : %12.3f +- %12.3f %s (%12.3f SIGMA)' % (error_value / measured_property.unit, error_uncertainty / measured_property.unit, str(measured_property.unit), relative_error)
            raise Exception(msg)
