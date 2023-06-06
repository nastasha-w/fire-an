import numpy as np
import numbers as num

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

def checkvalsmatch(val1, val2, verbose=False):
    if type(val1) != type(val2):
        if verbose:
            print(f'values {val1}, {val2} are of different types')
        return False
    elif isinstance(val1, num.Number):
        if isinstance(val1, int):
            same = val1 == val2
            if verbose and not same:
                print(f'Mismatch: {val1, val2}')
            return same
        else:
            same = np.isclose(val2, val2)
            if verbose and not same:
                print(f'Mismatch: {val1, val2}')
            return same
    elif isinstance(val1, type(np.string_(''))):
        same = val1.decode() == val2.decode()
        if verbose and not same:
            print(f'Mismatch: {val1, val2}')
        return same
    elif not hasattr(val1, '__len__'):
        same = val1 == val2
        if verbose and not same:
            print(f'Mismatch: {val1, val2}')
        return same
    else:
        if len(val1) != len(val2):
            if verbose:
                print(f'Array length mismatch:\n({len(val1)}) '
                      f'{val1}\n({len(val2)}) {val2}')
            return False
        else:
            if verbose:
                print(f'Recursing to check {val1}, {val2}')
            return np.all([checkvalsmatch(v1, v2, verbose=verbose) 
                           for v1, v2 in zip(val1, val2)])

def checkattrsmatch(grp1, grp2, verbose=False, subgroups=True):
    dct1 = {key: val for key, val in grp1.attrs.items()}
    dct2 = {key: val for key, val in grp2.attrs.items()}
    if set(dct1.keys()) != set(dct2.keys()):
        if verbose:
            msg = (f'for {grp1}, {grp2}, the sets of attributes were '
                   f'different:\n{set(dct1.keys())}\n{set(dct2.keys())}')
            print(msg)
        return False
    if subgroups:
        g1keys = set(grp1.keys())
        g2keys = set(grp2.keys())
        if not g1keys == g2keys:
            if verbose:
                msg = (f'for {grp1}, {grp2}, the sets of subgroups were '
                       f'different:\n{g1keys}\n{g2keys}')
                print(msg)
            return False
        else:
            gkeys = list(g1keys)
            subsame = np.all([checkattrsmatch(grp1[gkey], grp2[gkey], 
                                              verbose=verbose, subgroups=True)
                              for gkey in gkeys])
    else:
        subsame = True
    keys = list(dct1.keys())
    same = np.all([checkvalsmatch(dct1[key], dct2[key], verbose=verbose)
                   for key in keys])
    return same and subsame

def readgrp_todict(grp, subgroups=True):
    outdct = {}
    grpdct = {key: val for key, val in grp.attrs.items()}
    _keys = list(grpdct.keys())
    for key in _keys:
        if isinstance(grpdct[key], type(np.string_(''))):
            val = grpdct[key].decode()
            if val == 'None':
                outdct[key] = None
            elif val == 'dict':
                subkey = key + '_dict'
                outdct[key] = readgrp_todict(grp[subkey], 
                                             subgroups=subgroups)
            else:
                outdct[key] = val
        else:
            outdct[key] = grpdct[key]
    grpkeys = list(grp.keys())
    if subgroups:
        for key in grpkeys:
            outdct[key] = readgrp_todict(grp[key], subgroups=True)
    return outdct


