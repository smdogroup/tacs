"""
Sean Engelstad, January 2023
CAPS to TACS wrapper package
extends ESP/CAPS pyCAPS module to MPI for tacsAim
"""

# check the pyCAPS module can be loaded and only load remaining items if available
import importlib

caps_loader = importlib.util.find_spec("pyCAPS")

# main caps_struct module requires pyCAPS
if caps_loader is not None:
    from .tacs_model import *  # noqa: F403, E402

from .gif_writer import *  # noqa: F403, E402
from .analysis_function import *  # noqa: F403, E402
from .constraints import *  # noqa: F403, E402
from .egads_aim import *  # noqa: F403, E402
from .loads import *  # noqa: F403, E402
from .materials import *  # noqa: F403, E402
from .proc_decorator import *  # noqa: F403, E402
from .property import *  # noqa: F403, E402
from .tacs_aim import *  # noqa: F403, E402
from .variables import *  # noqa: F403, E402
