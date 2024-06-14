try: 
    from ._osfi_version import __version__, __author__
except ImportError:
    import warnings
    warnings.warn('No _osfi_version module found to be imported! Is OSFI running from the development tree?')
    __version__='(unknown)'
    __author__='(unknown)'

__all__ = ["CLP", "Logger", "Parameter", "ParamParserComplex", "ParamReader",
           "ParamType", "TimeValue", "UsageReader", "vt100", "__version__", "__author__"]
