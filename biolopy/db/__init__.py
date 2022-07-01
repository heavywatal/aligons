import os.path
from pathlib import Path

from ..util import cli

root = Path(os.path.expandvars(cli.config["db"]["root"])).expanduser()
