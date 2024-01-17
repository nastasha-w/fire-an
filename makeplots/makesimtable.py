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
               'm7e3': 7e3,
               'r7100': 7.1e3,
               'r3500': 3.5e3,
               'r4200': 4.2e3,
               'r28000': 3e4,
               'r57000': 6e4,
               }

# more sig. digits are listed for the FIRE-2 haloes
# but not generally for FIRE-3. Some other papers also give
# 1 sig. digit in their descriptions.
def get_resolution(simname):
    pt = sl.physlabel_from_simname(simname)
    if pt == 'FIRE-2':
        if simname.startswith('crheatfix'):
            respart = simname.split('_')[2]
        else:
            respart = simname.split('_')[1]
    elif pt == 'noBH-m12+':
        respart = simname.split('_')[3]
    else:
        respart = simname.split('_')[1]
    res = resolutions[respart]
    return res 

def get_refs(simname):
    return 'TODO'

sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2 \
          + sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
sims_f2md = sl.m12_f2md

def getdata(simname):
    _filldct = {}
    _filldct['ic'] = sl.ic_from_simname(simname)
    _filldct['phys'] = sl.plotlabel_from_physlabel[
        sl.physlabel_from_simname(simname)]
    res = get_resolution(simname)
    _filldct['gasres'] = res
    snap1 = max(sl.snaps_sr) if simname in sims_sr \
            else max(sl.snaps_hr) if simname in sims_hr \
            else max(sl.snaps_f2md) if simname in sims_f2md \
            else None
    snap0 = min(sl.snaps_sr) if simname in sims_sr \
            else min(sl.snaps_hr) if simname in sims_hr \
            else min(sl.snaps_f2md) if simname in sims_f2md \
            else None
    simpath = sl.dirpath_from_simname(simname)
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
    return _filldct

def maketable_main(simset='all'):
    if simset == 'FIRE-3':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-2':
        simnames = sl.m12_f2md
    elif simset == 'all':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'm12plus':
        simnames = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores

    simnames.sort()
    for sn in sl.buglist2:
        if sn in simnames:
            simnames.remove(sn)
    #ics = [sl.ic_from_simname(sn) for sn in simnames]
    #ics_clean = ['m12f', 'm13h113', 'm13h206'] 
    # ics_clean = [ic for ic in ics if sum([ic == _ic for _ic in ics]) == 3]
    #ics_clean = np.unique(ics_clean)
    #ics_clean.sort()
    #ics_rest = np.array(list(set(ics) - set(ics_clean)))
    #ics_rest.sort()
    #
    #simnames_clean = [simname for simname in simnames 
    #                  if sl.ic_from_simname(simname) in ics_clean]
    #simnames_clean.sort(key=lambda x: (sl.ic_from_simname(x), 
    #                                   sl.physlabel_from_simname(x)))
    #simnames_rest = [simname for simname in simnames 
    #                 if sl.ic_from_simname(simname) in ics_rest]
    #simnames_rest.sort(key=lambda x: (sl.ic_from_simname(x), 
    #                                  sl.physlabel_from_simname(x)))
    simgroupkeys = ['FIRE-2', 'noBH', 'noBH-m12+', 'AGN-noCR', 'AGN-CR']
    simgroups = {key: [sn for sn in simnames 
                       if sl.physlabel_from_simname(sn) == key]
                 for key in simgroupkeys}

    colsmain = ['{ic}', '{phys}', '{gasres}',
                '{mhalo0}', '{mstar0}', '{rhalo0}', 
                '{mhalo1}', '{mstar1}', '{rhalo1}']
    ncols = len(colsmain)
    aligndict = {'ic': 'l', 'phys': 'l', 'gasres': 'l',
                 'mhalo0': 'l', 'mstar0': 'l', 'rhalo0': 'l',
                 'mhalo1': 'l', 'mstar1': 'l', 'rhalo1': 'l',
                 'refs': 'c'}
    head0parts = ['', '', '',
                  '\\multicolumn{3}{c}{$z=1.0$}',
                  '\\multicolumn{3}{c}{$z=0.5$}',
                   ]
    head0 = ' \t & '.join(head0parts) + ' \t \\\\'
    head1dct = {'ic': 'ICs', 
                'phys': 'model', 
                'gasres': 'resolution',
                'mhalo0': '$\\mathrm{M}_{\\mathrm{vir}}$', 
                'mstar0': '$\\mathrm{M}_{\\star}$', 
                'rhalo0': '$\\mathrm{R}_{\\mathrm{vir}}$',
                'mhalo1': '$\\mathrm{M}_{\\mathrm{vir}}$', 
                'mstar1': '$\\mathrm{M}_{\\star}$', 
                'rhalo1': '$\\mathrm{R}_{\\mathrm{vir}}$',
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
    fmtdct = {'ic': '{ic}',
              'phys': '{phys}', 
              'gasres': '{gasres:.0e}',
              'mhalo0': '{mhalo0:.1e}',
              'mstar0': '{mstar0:.1e}',
              'rhalo0': '{rhalo0:.0f}',
              'mhalo1': '{mhalo1:.1e}',
              'mstar1': '{mstar1:.1e}',
              'rhalo1': '{rhalo1:.0f}',
              'refs': '{refs}',
              }
    hline = '\\hline'
    colspecfill = ' '.join(colsmain)
    colspec = colspecfill.format(**aligndict)
    start = f'\\begin{{tabular}}{{{colspec}}}'
    end = '\\end{tabular}'
    _fillmain = ' \t & '.join(colsmain) + ' \t \\\\'
    fillmain = _fillmain.format(**fmtdct)
    head1 = _fillmain.format(**head1dct)
    head2 = _fillmain.format(**head2dct)
    printlist = [start, hline, head0, head1, head2]
    
    for sgkey in simgroupkeys:
        simnames = simgroups[sgkey]
        if len(simnames) == 0:
            continue
        #shead = (f'\\multicolumn{{{ncols}}}{{c}}'
        #         f'{sl.plotlabel_from_physlabel(sgkey)} \\\\')
        printlist = printlist + [hline, hline]
        filldcts = {sn: getdata(sn) for sn in simnames}
        simnames.sort(key=lambda sn: filldcts[sn]['mhalo0'])
        for sn in simnames:
            _str = fillmain.format(**(filldcts[sn]))
            # 3e+05 -> 3e5, 3e+12 -> 3e12
            _str = _str.replace('+0', '')
            _str = _str.replace('+', '')
            printlist.append(_str)
        
    # for simname in simnames_clean:
    #     _filldct = {}
    #     _filldct['ic'] = sl.ic_from_simname(simname)
    #     _filldct['phys'] = sl.plotlabel_from_physlabel[
    #         sl.physlabel_from_simname(simname)]
    #     res = get_resolution(simname)
    #     _filldct['gasres'] = res
    #     snap1 = max(sl.snaps_sr) if simname in sims_sr \
    #             else max(sl.snaps_hr) if simname in sims_hr \
    #             else max(sl.snaps_f2md) if simname in sims_f2md \
    #             else None
    #     snap0 = min(sl.snaps_sr) if simname in sims_sr \
    #             else min(sl.snaps_hr) if simname in sims_hr \
    #             else min(sl.snaps_f2md) if simname in sims_f2md \
    #             else None
    #     simpath = sl.dirpath_from_simname(simname)
    #     _, _, halodat0 = cgp.readdata_cengalcen(simpath, snap0)
    #     _, _, halodat1 = cgp.readdata_cengalcen(simpath, snap1)
    #     _filldct['mhalo0'] = halodat0['halodata']['Mvir_g'] \
    #                          / c.solar_mass
    #     _filldct['rhalo0'] = halodat0['halodata']['Rvir_cm'] \
    #                          / (c.cm_per_mpc * 1e-3)
    #     _filldct['mstar0'] = halodat0['mstar_gal_g'] \
    #                          / c.solar_mass
    #     _filldct['mhalo1'] = halodat1['halodata']['Mvir_g'] \
    #                          / c.solar_mass
    #     _filldct['rhalo1'] = halodat1['halodata']['Rvir_cm'] \
    #                          / (c.cm_per_mpc * 1e-3)
    #     _filldct['mstar1'] = halodat1['mstar_gal_g'] \
    #                          / c.solar_mass
    #     _filldct['refs'] = get_refs(simname)
    #     _str = fillmain.format(**_filldct)
    #     # 3e+05 -> 3e5, 3e+12 -> 3e12
    #     _str = _str.replace('+0', '')
    #     _str = _str.replace('+', '')
    #     printlist.append(_str)
    # printlist = printlist + [hline, resthead, hline]
    # for simname in simnames_rest:
    #     _filldct = {}
    #     _filldct['ic'] = sl.ic_from_simname(simname)
    #     _filldct['phys'] = sl.plotlabel_from_physlabel[
    #         sl.physlabel_from_simname(simname)]
    #     res = get_resolution(simname)
    #     _filldct['gasres'] = res
    #     snap1 = max(sl.snaps_sr) if simname in sims_sr \
    #             else max(sl.snaps_hr) if simname in sims_hr \
    #             else max(sl.snaps_f2md) if simname in sims_f2md \
    #             else None
    #     snap0 = min(sl.snaps_sr) if simname in sims_sr \
    #             else min(sl.snaps_hr) if simname in sims_hr \
    #             else min(sl.snaps_f2md) if simname in sims_f2md \
    #             else None
    #     simpath = sl.dirpath_from_simname(simname)
    #     _, _, halodat0 = cgp.readdata_cengalcen(simpath, snap0)
    #     _, _, halodat1 = cgp.readdata_cengalcen(simpath, snap1)
    #     _filldct['mhalo0'] = halodat0['halodata']['Mvir_g'] \
    #                          / c.solar_mass
    #     _filldct['rhalo0'] = halodat0['halodata']['Rvir_cm'] \
    #                          / (c.cm_per_mpc * 1e-3)
    #     _filldct['mstar0'] = halodat0['mstar_gal_g'] \
    #                          / c.solar_mass
    #     _filldct['mhalo1'] = halodat1['halodata']['Mvir_g'] \
    #                          / c.solar_mass
    #     _filldct['rhalo1'] = halodat1['halodata']['Rvir_cm'] \
    #                          / (c.cm_per_mpc * 1e-3)
    #     _filldct['mstar1'] = halodat1['mstar_gal_g'] \
    #                          / c.solar_mass
    #     _filldct['refs'] = get_refs(simname)
    #     _str = fillmain.format(**_filldct)
    #     # 3e+05 -> 3e5, 3e+12 -> 3e12
    #     _str = _str.replace('+0', '')
    #     _str = _str.replace('+', '')
    #     printlist.append(_str)
    printlist = printlist + [hline, end]
    table = '\n'.join(printlist)
    print(table)

def getlongname(simname):
    if sl.physlabel_from_simname(simname) == 'FIRE-2':
        if simname.startswith('crheatfix'):
            longname = '_'.join(simname.split('_')[1:])
            longname = 'core + metal diffusion + CR heating fix, ' \
                       + longname
        else:
            longname = 'core + metal diffusion, ' + simname
    elif sl.physlabel_from_simname(simname) == 'noBH-m12+':
        if simname.startswith('fire3nobh_plus_'):
            longname = '_'.join(simname.split('_')[2:])
            longname = 'high-mass m12 FIRE-3 NoBH, ' + longname
    else:
        longname = simname
    longname = longname.replace('_', '\\_')
    return longname

def maketable_appendix(simset='all'):
    if simset == 'FIRE-3':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-2':
        simnames = sl.m12_f2md
    elif simset == 'all':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md + \
                   sl.m13_hr_all2 + sl.m13_sr_all2

    simnames.sort()
    for sn in sl.buglist2:
        if sn in simnames:
            simnames.remove(sn)
    ics = [sl.ic_from_simname(sn) for sn in simnames]
    ics_clean = ['m12f', 'm13h113', 'm13h206'] 
    # ics_clean = [ic for ic in ics if sum([ic == _ic for _ic in ics]) == 3]
    ics_clean = np.unique(ics_clean)
    ics_clean.sort()
    ics_rest = np.array(list(set(ics) - set(ics_clean)))
    ics_rest.sort()
    
    simnames_clean = [simname for simname in simnames 
                      if sl.ic_from_simname(simname) in ics_clean]
    simnames_clean.sort(key=lambda x: (sl.ic_from_simname(x), 
                                       sl.physlabel_from_simname(x)))
    simnames_rest = [simname for simname in simnames 
                     if sl.ic_from_simname(simname) in ics_rest]
    simnames_rest.sort(key=lambda x: (sl.ic_from_simname(x), 
                                      sl.physlabel_from_simname(x)))
    colsmain = ['{ic}', '{phys}', '{simname}']
    ncols = len(colsmain)
    aligndict = {'ic': 'l', 'phys': 'l', 'simname': 'l'}
    head1dct = {'ic': 'ICs', 
                'phys': 'model', 
                'simname': 'FIRE collaboration internal name',
                }
    fmtdct = {'ic': '{ic}',
              'phys': '{phys}', 
              'simname': '{simname}'
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
    printlist = [start, hline, head1, hline, cleanhead, hline]

    for simname in simnames_clean:
        _filldct = {}
        _filldct['ic'] = sl.ic_from_simname(simname)
        _filldct['phys'] = sl.plotlabel_from_physlabel[
            sl.physlabel_from_simname(simname)]
        _filldct['simname'] = getlongname(simname)
        _str = fillmain.format(**_filldct)
        printlist.append(_str)
    printlist = printlist + [hline, resthead, hline]
    for simname in simnames_rest:
        _filldct = {}
        _filldct['ic'] = sl.ic_from_simname(simname)
        _filldct['phys'] = sl.plotlabel_from_physlabel[
            sl.physlabel_from_simname(simname)]
        _filldct['simname'] = getlongname(simname)
        _str = fillmain.format(**_filldct)
        printlist.append(_str)
    printlist = printlist + [hline, end]
    table = '\n'.join(printlist)
    print(table)

def printsims_used(simset='all'):
    if simset == 'FIRE-3':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-2':
        simnames = sl.m12_f2md
    elif simset == 'all':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-3 m12+':
        simnames = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    
    simnames.sort()
    for sn in sl.buglist2:
        if sn in simnames:
            simnames.remove(sn)
    simnames.sort(key=lambda x: (sl.ic_from_simname(x), 
                                 sl.physlabel_from_simname(x)))
    for simname in simnames:
        _sn = getlongname(simname)
        _sn = _sn.replace('\_', '_') # remove LaTeX escapes
        print(_sn)


def print_checklist(simset='FIRE-3'):
    if simset == 'FIRE-3':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-2':
        simnames = sl.m12_f2md
    elif simset == 'all':
        simnames = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md + \
                   sl.m13_hr_all2 + sl.m13_sr_all2
    elif simset == 'FIRE-3 m12+':
        simnames = sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    sims_md = sl.m12_f2md
    sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2 \
              + sl.m12plus_f3nobh + sl.m12plus_f3nobh_lores
    sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2
    for simname in simnames:
        _path = '.' + sl.dirpath_from_simname(simname) + 'output/'
        if simname in sims_md:
            snapshots = sl.snaps_f2md
        elif simname in sims_hr:
            snapshots = sl.snaps_hr
        elif simname in sims_sr:
            snapshots = sl.snaps_sr

        print(_path + ' \t ' + str(snapshots))