

import pandas as pd
# Things that are specific to python 2.x
import sys
if sys.version.startswith('2'):
    import pickle
    # custom magics for pasting
elif sys.version.startswith('3'):
    import importlib
    from importlib import reload