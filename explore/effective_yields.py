'''
Calculate total stellar masses and metal masses in
zoom simulation volumes
'''

import h5py
import matplotlib.gridspec as gsp
import matplotlib.pyplot as plt
from matplotlib.transforms import Transform
import numpy as np

import fire_an.mainfunc.get_qty as gq
import fire_an.readfire.readin_fire_data as rfd
import fire_an.simlists as sl
import fire_an.utils.h5utils as h5u

snaps_hr = sl.snaps_hr
snaps_sr = sl.snaps_sr
snaps_md = sl.snaps_f2md
sims_hr = sl.m12_hr_all2 + sl.m13_hr_all2
sims_sr = sl.m12_sr_all2 + sl.m13_sr_all2 + sl.m12_fire3x_tests
sims_md = sl.m12_f2md
snaps_z = sl.snaps_z
physcolors = sl.physcolors.copy()
physcolors.update({'FIRE-3x-scmodules': sl._physcolors.cyan,
                   'FIRE-3x-constpterm': sl._physcolors.purple})


def get_totals(simname, snapnum, outname):
    metals = ['total', 'Oxygen', 'Neon', 'Carbon', 'Nitrogen', 'Iron',
              'Magnesium', 'Sulfur']
    maptypes = ['Mass'] + ['Metal'] * len(metals)
    maptype_argss = [{}] + [{'element': val, 'density': False}
                            for val in metals]

    dirpath = sl.dirpath_from_simname(simname)
    snap = rfd.get_Firesnap(dirpath, snapnum)
    with h5py.File(outname, 'a') as f:
        hed = f.create_group('Header')
        hed.attrs.create('simname', np.string_(simname))
        hed.attrs.create('dirpath', np.string_(dirpath))
        hed.attrs.create('snapnum', snapnum)
        cosmopars = snap.cosmopars.getdct()
        csm = hed.create_group('cosmopars')
        h5u.savedict_hdf5(csm, cosmopars)
        f.create_group('gas')
        f.create_group('stars')

        for parttype, ptlabel in [(0, 'gas'), (4, 'stars')]:
            for maptype, maptype_args in zip(maptypes, maptype_argss):
                qty, toCGS, todoc = gq.get_qty(snap, parttype, 
                                               maptype, maptype_args,
                                               filterdct=None)
                tot = np.sum(qty)

                grpname = ('Mass' if maptype == 'Mass' 
                           else 'MetalMass_' + maptype_args['element'])
                cur = f[ptlabel].create_group(grpname)
                cur.attrs.create('total mass', tot)
                cur.attrs.create('toCGS', toCGS)
                cdoc = cur.create_group('doc')
                h5u.savedict_hdf5(cdoc, todoc)

def run_totals(index):
    outdir = '/projects/b1026/nastasha/hists/gas_stars_metals/'
    # leaving out the fire3_m12plus halos for npw
    if index >= 0 and index < 108:
        ind = index - 0
        simnames = sl.m12_hr_all2 # 18
        snapnums = sl.snaps_hr # 6
    elif index >= 108 and index < 132:
        ind = index - 108
        simnames = sl.m12_sr_all2 # 4
        snapnums = sl.snaps_sr # 6 
    elif index >= 132 and index < 144:
        ind = index - 132
        simnames = sl.m13_hr_all2 # 2
        snapnums = sl.snaps_hr # 6
    elif index >= 144 and index < 234:
        ind = index - 144
        simnames = sl.m13_sr_all2 # 15
        snapnums = sl.snaps_sr # 6
    elif index >= 234 and index < 294:
        ind = index - 234
        simnames = sl.m12_f2md # 10
        snapnums = sl.snaps_f2md #6
    elif index >= 294 and index < 318:
        # frontera output dir
        outdir = '/scratch1/08466/tg877653/output/hists/gas_stars_metals/'
        ind = index - 294
        simnames = sl.m12_fire3x_tests # 4
        snapnums = sl.snaps_sr #6
    
    nmi = ind // len(snapnums)
    sni = ind % len(snapnums)
    simname = simnames[nmi]
    snapnum = snapnums[sni]
    
    outname = outdir + f'total_MassZ_stars_gas_{simname}_{snapnum}.hdf5'

    get_totals(simname, snapnum, outname)
            
def get_yielddata(simname, species='total', source='all'):
    snapnums = snaps_hr if simname in sims_hr \
               else snaps_sr if simname in sims_sr \
               else snaps_md if simname in sims_md \
               else None
    zs = []
    yields = []
    if source == 'all':
        sources = ['gas', 'stars']
    else:
        sources = [source]
    for i, snapnum in enumerate(snapnums):
        filen = ('/projects/b1026/nastasha/hists/gas_stars_metals/'
                 f'total_MassZ_stars_gas_{simname}_{snapnum}.hdf5')
        zkey = snaps_z[i]
        zs.append(zkey)
        with h5py.File(filen, 'r') as f:
            mstar = f['stars/Mass'].attrs('total mass')
            conv = f['stars/Mass'].attrs('toCGS')
            mstar *= conv
            mZ = 0.
            for source in sources:
                path = f'{source}/MetalMass_{species}'
                _mZ = f[path].attrs('total mass')
                _conv = f[path].attrs('toCGS')
                if species != 'total':
                    _conv *= gq.elt_atomw_cgs(species)
                mZ += _mZ * _conv
        yields.append(mZ / mstar)
    return zs, yields
    

def plotyields(species='total', source='all'):
    simnames_all = sl.m12_hr_all2 + sl.m13_hr_all2 \
                   + sl.m12_sr_all2 + sl.m13_sr_all2 \
                   + sl.m12_fire3x_tests + sl.m12_f2md
    bugsims = sl.buglist2
    ics = [sl.ic_from_simname(sn) for sn in simnames_all]
    ics_all = np.unique(ics)
    ics_all.sort()
    ics_special = ['m12f', 'm12m', 'm12i']
    allyields_bug = {}
    specialyields_bug = {}
    allyields_nobug = {}
    specialyields_nobug = {}
    
    fontsize = 12
    npanels = len(ics_all)
    ncols = min(npanels, 4)
    nrows = (npanels - 1) // ncols + 1
    _nrows = nrows + 1
    panelsize = 2.5
    hspace = 0.4
    wspace = 0.
    width = panelsize * ncols * (1. + (ncols - 1.) * wspace / ncols)
    height = panelsize * _nrows * (1. + hspace / _nrows)

    fig = plt.figure(figsize=(width, height))
    _grid = gsp.GridSpec(nrows=2, ncols=1, 
                         height_ratios=[nrows, 1], hspace=hspace)
    grid = gsp.GridSpecFromSubplotSpec(nrows=nrows, ncols=ncols, 
                                       hspace=0., wspace=0.,
                                       subplot_spec=_grid[0])
    mainaxes = [fig.add_subplot(grid[i // ncols, i % ncols]) 
                for i in range(npanels)]
    hgrid = gsp.GridSpecFromSubplotSpec(nrows=1, ncols=3, 
                                        hspace=0., wspace=0.,
                                        subplot_spec=_grid[1])
    haxes = [fig.add_subplot(hgrid[0, i]) 
             for i in range(3)]
    
    title1 = 'total metal mass' if species =='total' \
             else f'{species} mass'
    title2 = ' in stars and gas' if source == 'all' \
             else f' in {source}'
    title = title1 + title2 + ' / current stellar mass'
    fig.suptitle(title)
    
    for axi, ic in enumerate(ics_all):
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        ax = mainaxes[axi]
        ax.tick_params(which='both', labelsize=fontsize - 1,
                       direction='in', right=True, top=True,
                       grid=True, labelbottom=dobottom,
                       labelleft=doleft)
        if dobottom:
            ax.set_xlabel('redshift', fontsize=fontsize)
        if doleft:
            ax.set_ylabel('$\\mathrm{M}_{\\mathrm{Z}} \\,/\\,'
                          '\\mathrm{M}_{\\star}$')
        ax.text(0.95, 0.95, ic, fontsize=fontsize - 1.,
                transform=ax.transAxes,
                horizontalalignment='right',
                verticalalignment='top')
        for simname in simnames_all:
            if sl.ic_from_simname(simname) != ic:
                continue
            phys = sl.physlabel_from_simname(simname)
            hasbug = simname in bugsims
            color = physcolors[phys]
            linestyle = 'dashed' if hasbug else 'solid'
            zs, yields = get_yielddata(simname, species=species,
                                       source=source)
            ax.plot(zs, yields, color=color, linestyle=linestyle,
                    linewidth=1.5, marker='o')
            if hasbug:
                if ic in ics_special:
                    if phys in specialyields_bug:
                        specialyields_bug[phys] += yields
                    else:
                        specialyields_bug[phys] = yields
                if phys in allyields_bug:
                    allyields_bug[phys] += yields
                else:
                    allyields_bug[phys] = yields
            else:
                if ic in ics_special:
                    if phys in specialyields_nobug:
                        specialyields_nobug[phys] += yields
                    else:
                        specialyields_nobug[phys] = yields
                if phys in allyields_nobug:
                    allyields_nobug[phys] += yields
                else:
                    allyields_nobug[phys] = yields
                




    