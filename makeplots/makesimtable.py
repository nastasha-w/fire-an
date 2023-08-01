import numpy as np

import fire_an.mainfunc.cengalprop as cgp
import fire_an.simlists as sl
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol

# first go; check actual values, Msun or Msun/h
resolutions = {'m3e5': 3e5,
               'm3e4': 3e4,
               'm6e4': 6e4,
               'm4e3': 4e3,
               'm7e3': 7e3}

def get_model(simname):
    if '_sdp1e10_' in simname:
        return 'noBH'
    elif '_MHDCRspec1_' in simname:
        return 'AGN-CR'
    else:
        return 'AGN-noCR'

def get_ic(simname):
    return simname.split('_')[0]

# more sig. digits are listed for the FIRE-2 haloes
# but not generally for FIRE-3. Some other papers also give
# 1 sig. digit in their descriptions.
def get_resolution(simname):
    respart = simname.split('_')[1]
    res = resolutions[respart]
    return res 

def get_refs(simname):
    return 'TODO'

def getpath(simname):
    dp2 = '_'.join(simname.split('_')[:2])
    if dp2.startswith('m13h02_'):
        dp2 = dp2.replace('m13h02', 'm13h002')
    simpath = '/'.join([ol.simdir_fire, dp2, simname]) 
    return simpath

def maketable_main(simset='FIRE-3'):
    if simset == 'FIRE-3':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + \
                sl.m13_hr_all2 + sl.m13_sr_all2
        sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
        sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-2':
        simnames = sl.m12_f2md
    simnames.sort()
    for sn in sl.buglist1:
        if sn in simnames:
            simnames.remove(sn)
    ics = [get_ic(sn) for sn in simnames]
    ics_clean = ['m12f', 'm13h113', 'm13h206'] 
    # ics_clean = [ic for ic in ics if sum([ic == _ic for _ic in ics]) == 3]
    ics_clean = np.unique(ics_clean)
    ics_clean.sort()
    ics_rest = np.array(list(set(ics) - set(ics_clean)))
    ics_rest.sort()
    
    simnames_clean = [simname for simname in simnames 
                      if get_ic(simname) in ics_clean]
    simnames_clean.sort(key=lambda x: (get_ic(x), get_model(x)))
    simnames_rest = [simname for simname in simnames 
                     if get_ic(simname) in ics_rest]
    simnames_rest.sort(key=lambda x: (get_ic(x), get_model(x)))

    colsmain = ['{ic}', '{phys}', '{gasres}',
                '{mhalo0}', '{mstar0}', '{rhalo0}', 
                '{mhalo1}', '{mstar1}', '{rhalo1}']
    ncols = len(colsmain)
    aligndict = {'ic': 'l', 'phys': 'l', 'gasres': 'r',
                 'mhalo0': 'r', 'mstar0': 'r', 'rhalo0': 'r',
                 'mhalo1': 'r', 'mstar1': 'r', 'rhalo1': 'r',
                 'refs': 'c'}
    head1dct = {'ic': 'ICs', 'phys': 'model', 'gasres': 'resolution',
                'mhalo0': '$\\mathrm{M}_{\\mathrm{vir}}(1.0)$', 
                'mstar0': '$\\mathrm{M}_{\\star}(1.0)$', 
                'rhalo0': '$\\mathrm{R}_{\\mathrm{vir}}(1.0)$',
                'mhalo1': '$\\mathrm{M}_{\\mathrm{vir}}(0.5)$', 
                'mstar1': '$\\mathrm{M}_{\\star}(0.5)$', 
                'rhalo1': '$\\mathrm{R}_{\\mathrm{vir}}(0.5)$',
                'refs': 'references',
                }
    head2dct = {'ic': '', 'phys': '', 
                'gasres': ('$[\\mathrm{M}_{\\odot}]$'),
                'mhalo0': ('$[\\mathrm{M}_{\\odot}]$'),
                'mstar0': ('$[\\mathrm{M}_{\\odot}]$'),
                'rhalo0': '[pkpc]',
                'mhalo1': ('$[\\mathrm{M}_{\\odot}]$'),
                'mstar1': ('$[\\mathrm{M}_{\\odot}]$'),
                'rhalo1': '[pkpc]',
                'refs': '',
                }
    fmtdct = {'ic': '{ic}', 'phys': '{phys}', 
              'gasres': '{gasres:.1e}',
              'mhalo0': '{mhalo0:.1e}',
              'mstar0': '{mstar0:.1e}',
              'rhalo0': '{rhalo0:.0f}',
              'mhalo1': '{mhalo1:.1e}',
              'mstar1': '{mstar1:.1e}',
              'rhalo1': '{rhalo1:.0f}',
              'refs': '{refs}',
              }
    cleanhead = (f'\\multicolumn{{{ncols}}}{{c}}'
                 '{clean sample} \\\\')
    resthead = (f'\\multicolumn{{{ncols}}}{{c}}'
                 '{full sample} \\\\')
    hline = '\\hline'
    colspecfill = ' '.join(colsmain)
    colspec = colspecfill.format(**aligndict)
    start = f'\\begin{{tabular}}{{{colspec}}}'
    end = '\\end{tabular}'
    _fillmain = ' \t & '.join(colsmain) + ' \t \\\\'
    fillmain = _fillmain.format(**fmtdct)
    head1 = _fillmain.format(**head1dct)
    head2 = _fillmain.format(**head2dct)
    printlist = [start, hline, head1, head2, hline, cleanhead, hline]

    for simname in simnames_clean:
        _filldct = {}
        _filldct['ic'] = get_ic(simname)
        _filldct['phys'] = get_model(simname)
        res = get_resolution(simname)
        _filldct['gasres'] = res
        snap1 = max(sl.snaps_sr) if simname in sims_sr \
                else max(sl.snaps_hr) if simname in sims_hr \
                else None
        snap0 = min(sl.snaps_sr) if simname in sims_sr \
                else min(sl.snaps_hr) if simname in sims_hr \
                else None
        simpath = getpath(simname)
        _, _, halodat0 = cgp.readdata_cengalcen(simpath, snap0)
        _, _, halodat1 = cgp.readdata_cengalcen(simpath, snap1)
        _filldct['mhalo0'] = halodat0['halodata']['Mvir_g'] \
                             / c.solar_mass
        _filldct['rhalo0'] = halodat0['halodata']['Rvir_cm'] \
                             / (c.cm_per_mpc * 1e-3)
        _filldct['mstar0'] = halodat0['mstar_gal_g'] \
                             / c.solar_mass
        _filldct['mhalo1'] = halodat1['halodata']['Mvir_g'] \
                             / c.solar_mass
        _filldct['rhalo1'] = halodat1['halodata']['Rvir_cm'] \
                             / (c.cm_per_mpc * 1e-3)
        _filldct['mstar1'] = halodat1['mstar_gal_g'] \
                             / c.solar_mass
        _filldct['refs'] = get_refs(simname)
        printlist.append(fillmain.format(**_filldct))
    printlist = printlist + [hline, resthead, hline]
    for simname in simnames_rest:
        _filldct = {}
        _filldct['ic'] = get_ic(simname)
        _filldct['phys'] = get_model(simname)
        res = get_resolution(simname)
        _filldct['gasres'] = res
        snap1 = max(sl.snaps_sr) if simname in sims_sr \
                else max(sl.snaps_hr) if simname in sims_hr \
                else None
        snap0 = min(sl.snaps_sr) if simname in sims_sr \
                else min(sl.snaps_hr) if simname in sims_hr \
                else None
        simpath = getpath(simname)
        _, _, halodat0 = cgp.readdata_cengalcen(simpath, snap0)
        _, _, halodat1 = cgp.readdata_cengalcen(simpath, snap1)
        _filldct['mhalo0'] = halodat0['halodata']['Mvir_g'] \
                             / c.solar_mass
        _filldct['rhalo0'] = halodat0['halodata']['Rvir_cm'] \
                             / (c.cm_per_mpc * 1e-3)
        _filldct['mstar0'] = halodat0['mstar_gal_g'] \
                             / c.solar_mass
        _filldct['mhalo1'] = halodat1['halodata']['Mvir_g'] \
                             / c.solar_mass
        _filldct['rhalo1'] = halodat1['halodata']['Rvir_cm'] \
                             / (c.cm_per_mpc * 1e-3)
        _filldct['mstar1'] = halodat1['mstar_gal_g'] \
                             / c.solar_mass
        _filldct['refs'] = get_refs(simname)
        printlist.append(fillmain.format(**_filldct))
    printlist = printlist + [hline, end]
    table = '\n'.join(printlist)
    print(table)
        


