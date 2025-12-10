import importlib.metadata

assert __package__
__doc__ = importlib.metadata.metadata(__package__)["Summary"]
__version__ = importlib.metadata.version(__package__)
