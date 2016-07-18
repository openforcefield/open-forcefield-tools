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

#=============================================================================================
# DATASET
#=============================================================================================

class PhysicalPropertyDataset(object):
    """A dataset of physical property measurements.

    TODO:
    * Should this subclass list or set?
    """
    pass

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
