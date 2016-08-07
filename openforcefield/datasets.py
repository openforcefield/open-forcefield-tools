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
* Create PhysicalPropertyDataset API.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy

#=============================================================================================
# DATASET
#=============================================================================================

class PhysicalPropertyDataset(object):
    """A dataset of physical property measurements.

    Implements the container API.

    Properties
    ----------

    """
    def __init__(self, dataset=None):
        """Create a physical property dataset.

        Parameters
        ----------
        dataset : iterable, optional, default=None
            If a dataset is specified, it is copied to create a new dataset.

        """
        self._measurements = list()

        if dataset is not None:
            for measurement in dataset:
                self._measurements.append( copy.deepcopy(measurement) )

    def __len__(self):
        return len(self._measurements)

    def __getitem__(self, doi):
        """Return all measurements with the given DOI.
        """
        measurements = list()
        for measurement in self._measurements:
            if measurement.doi == doi:
                measurements.append(measurement)
        if len(measurements) == 0:
            return KeyError("DOI (%s) not found in dataset." % doi)
        return measurements

    def __delitem__(self, doi):
        """Delete all items with the given DOI.
        """
        measurement = self[doi]
        self._measurements.remove(measurement)

    def __iter__(self):
        for measurement in self._measurements:
            yield measurement

    def __contains__(self, measurement):
        for measurement2 in self._measurements:
            if measurement == measurement2:
                return True
        return False

#=============================================================================================
# THERMOML DATASET
#=============================================================================================

class ThermoMLDataset(PhysicalPropertyDataset):
    """A dataset of physical property measurements created from a ThermoML dataset.

    Examples
    --------

    For example, we can use the DOI `10.1016/j.jct.2005.03.012` as a key for retrieving the dataset from the ThermoML Archive:

    >>> dataset = ThermoMLDataset('10.1016/j.jct.2005.03.012')

    You can also specify multiple ThermoML Archive keys to create a dataset from multiple ThermoML files:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataset(thermoml_keys)

    You can see which DOIs contribute to the current `ThermoMLDataset` with the convenience functions:

    >>> thermoml_keys = ['10.1021/acs.jced.5b00365', '10.1021/acs.jced.5b00474']
    >>> dataset = ThermoMLDataset(thermoml_keys)

    """
    def __init__(self, doi=None, url=None):
        """Retrieve ThermoML files matching specified keys from the specified URL.

        Parameters
        ----------
        doi : str or list of str, optional, default=None
            If specified, ThermoML files with the specified DOI keys will be retrieved
        url : str, optional, default=None
            If specified, this URL (which may point to a local filesystem) will be used instead of the ThermoML Archive.

        """
        pass
