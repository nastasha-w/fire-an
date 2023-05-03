import numpy as np

def savedict_hdf5(grp, dct):
    for key in dct:
        val = dct[key]
        if isinstance(val, type('')):
            val = np.string_(val)
        elif val is None:
            val = np.string_('None')
        elif isinstance(val, dict):
            sgrp = grp.create_group(key + '_dict')
            _val = val.copy()
            savedict_hdf5(sgrp, _val)
            val = np.string_('dict')
        grp.attrs.create(key, val)