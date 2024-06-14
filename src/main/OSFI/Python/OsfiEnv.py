#
# openSF Integration Libraries (OSFI)
# Deimos Space, S.L.U.
# 
# This file is part of OSFI. OSFI is free software; you can redistribute it
# and/or modify it under the terms of the 'ESA Software Community Licence Permissive' as
# published by the European Space Agency; either version 2.4 of the License,
# or (at your option) any later version. You should have received a
# copy of the 'ESA Software Community Licence Permissive - v2.4' along with this program
# or one can be found at <http://eop-cfi.esa.int/index.php/docs-and-mission-data/licensing-documents>.
#

class OsfiEnv:
    """
    Singleton store of environmental settings for OSFI.
    Created so that init-once and environment-based configurations are easier to mock
    for unit-testing purposes.
    """

    def __init__(self, base_cfg=None):
        self._base_cfg = base_cfg # Null for normal operation, redirecting external_cfg to os.environ
        self._values = {}

    def external_var(self, key):
        """
        Retrieve a value from the external configuration, which is normally os.environ.
        :param key: Name of the environment variable.
        :returns: Value, or None if the variable is not present.
        """
        if self._base_cfg is not None:
            return self._base_cfg.get(key, None)
        else:
            import os
            return os.environ.get(key, None)

    # Methods delegated to the _values dict
    def __getitem__(self, key):
        """Retrieve value from the shared storage, or None if no value is present under that key."""
        return self._values.get(key, None)
    def __setitem__(self, key, value):
        """Store a value under a given ID"""
        self._values[key] = value
    def __delitem__(self, key):
        del self._values[key]
    def __contains__(self, key):
        return key in self._values

# Shared instance. In normal operation this is only written to once, but it can be replaced by UTs
OsfiEnv._INSTANCE = OsfiEnv()
def get_env():
    """Retrieve the shared singleton instance of the OSFI environment."""
    return OsfiEnv._INSTANCE
