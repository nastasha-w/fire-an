'''
Calculate total stellar masses and metal masses in
zoom simulation volumes
'''

import h5py
import matplotlib.gridspec as gsp
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
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
            mstar = f['stars/Mass'].attrs['total mass']
            conv = f['stars/Mass'].attrs['toCGS']
            mstar *= conv
            mZ = 0.
            for source in sources:
                path = f'{source}/MetalMass_{species}'
                _mZ = f[path].attrs['total mass']
                _conv = f[path].attrs['toCGS']
                if species != 'total':
                    _conv *= gq.elt_atomw_cgs(species)
                mZ += _mZ * _conv
        yields.append(mZ / mstar)
    return zs, yields
    

def plotyields(species='total', source='all', massset='m12',
               outname=None):
    simnames_all = sl.m12_hr_all2 + sl.m13_hr_all2 \
                   + sl.m12_sr_all2 + sl.m13_sr_all2 \
                   + sl.m12_fire3x_tests + sl.m12_f2md
    bugsims = sl.buglist2
    ics = [sl.ic_from_simname(sn) for sn in simnames_all
           if sn.startswith(massset)]
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
    hspace = 0.2
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
                                        hspace=0., wspace=0.2,
                                        subplot_spec=_grid[1])
    haxes = [fig.add_subplot(hgrid[0, i]) 
             for i in range(3)]
    
    title1 = 'total metal mass' if species =='total' \
             else f'{species} mass'
    title2 = ' in stars and gas' if source == 'all' \
             else f' in {source}'
    title = title1 + title2 + ' / current stellar mass'
    fig.suptitle(title)

    yieldlabel = '$\\mathrm{M}_{\\mathrm{Z}} \\,/\\, \\mathrm{M}_{\\star}$'
    
    for axi, ic in enumerate(ics_all):
        dobottom = axi >= npanels - ncols
        doleft = axi % ncols == 0
        ax = mainaxes[axi]
        ax.tick_params(which='both', labelsize=fontsize - 1,
                       direction='in', right=True, top=True,
                       labelbottom=dobottom,
                       labelleft=doleft)
        ax.grid(visible=True, which='major', axis='both')
        if dobottom:
            ax.set_xlabel('redshift', fontsize=fontsize)
        if doleft:
            ax.set_ylabel(yieldlabel, fontsize=fontsize)
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
    ylims = [ax.get_ylim() for ax in mainaxes]
    ymin = min([l[0] for l in ylims])
    ymax = max([l[1] for l in ylims])
    [ax.set_ylim(ymin, ymax) for ax in mainaxes]
    
    yieldrange = ymax - ymin
    nbins = 20
    binsize = yieldrange / nbins
    physkeys = set(allyields_bug.keys()) | set(allyields_nobug.keys()) \
               | set(specialyields_bug.keys()) \
               | set(specialyields_nobug.keys())
    physkeys = list(physkeys)
    nkeys = len(physkeys)
    offset = 0.5 * binsize / nkeys
    bins = {key: np.linspace(ymin - i * offset,
                             ymax + binsize - i * offset,
                             nbins + 1)
            for i, key in enumerate(physkeys)}
    hax_all = haxes[0]
    hax_all.tick_params(which='both', labelsize=fontsize - 1,
                        direction='in', right=True, top=True,
                        labelbottom=True, labelleft=True)
    hax_spec = haxes[1]
    hax_spec.tick_params(which='both', labelsize=fontsize - 1,
                        direction='in', right=True, top=True,
                        labelbottom=True, labelleft=True)
    hax_all.set_xlabel(yieldlabel, fontsize=fontsize)
    hax_spec.set_xlabel(yieldlabel, fontsize=fontsize)
    hax_all.set_ylabel('snapshot, halo count', fontsize=fontsize)
    hax_all.text(0.95, 0.95, 'all halos',
                 fontsize=fontsize - 1,
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=hax_all.transAxes)
    hax_spec.text(0.95, 0.95, ', '.join(ics_special),
                  fontsize=fontsize - 1,
                  horizontalalignment='right',
                  verticalalignment='top',
                  transform=hax_spec.transAxes)
    for key in physkeys:
        if key in allyields_nobug:
            hax_all.hist(allyields_nobug[key], bins=bins[key],
                        color=physcolors[key], linestyle='solid',
                        alpha=0.7, histtype='step', linewidth=1.2)
        if key in allyields_bug:
            hax_all.hist(allyields_bug[key], bins=bins[key],
                        color=physcolors[key], linestyle='dashed',
                        alpha=1., histtype='step', linewidth=1.5)
        if key in specialyields_nobug:
            hax_spec.hist(specialyields_nobug[key], bins=bins[key],
                        color=physcolors[key], linestyle='solid',
                        alpha=0.7, histtype='step', linewidth=1.2)
        if key in specialyields_bug:
            hax_spec.hist(specialyields_bug[key], bins=bins[key],
                        color=physcolors[key], linestyle='dashed',
                        alpha=1., histtype='step', linewidth=1.5)
    
    lax = haxes[2]
    lax.axis('off')
    physkeys.sort()
    leglines = [mlines.Line2D((), (), linestyle='solid', linewidth=1.5,
                              color=physcolors[key], 
                              label=sl.plotlabel_from_physlabel[key])
                for key in physkeys]
    leglines = leglines + \
               [mlines.Line2D((), (), linestyle='solid', linewidth=1.5,
                              color='gray', 
                              label='no bug'),
                mlines.Line2D((), (), linestyle='dashed', linewidth=1.5,
                              color='gray', 
                              label='bug'), 
                ]
    lax.legend(handles=leglines, fontsize=fontsize - 1)
    
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def runplots_effyields():
    outdir = '/projects/b1026/nastasha/imgs/effective_yields/'

    for species in ['total', 'Oxygen', 'Neon', 'Carbon', 'Nitrogen', 
                    'Iron', 'Magnesium', 'Sulfur']:
        for source in ['all', 'stars', 'gas']:
            for massset in ['m12', 'm13']:
                outname = (f'effective_yields_{massset}_{species}'
                           f'_in_{source}.pdf')
                plotyields(species=species, source=source, massset=massset,
                           outname=outdir + outname)
        
                




    