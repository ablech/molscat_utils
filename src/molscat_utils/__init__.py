# Expose submodules for easy import
from . import (
    wrapper,
    parser,
    cli,
)

__version__ = "0.1.0.dev0"

import git
try:
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha[0:7]
except:
    sha = "unknown"
__git_version__ = __version__ + f"+{sha}"