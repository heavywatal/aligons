import os.path
from pathlib import Path

from ..util import config

root = Path(os.path.expandvars(config["db"]["root"])).expanduser()
