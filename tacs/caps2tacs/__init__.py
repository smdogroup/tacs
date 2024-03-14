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
    from .tacs_model import *

openmdao_loader = importlib.util.find_spec("openmdao")

# import the openmdao component only if available to load
if openmdao_loader is not None:
    from .tacs_component import *

from .gif_writer import *
from .analysis_function import *
from .constraints import *
from .egads_aim import *
from .loads import *
from .materials import *
from .proc_decorator import *
from .property import *
from .tacs_aim import *
from .variables import *
