
import h5py
import numpy as np

import matplotlib.collections as mcol
import matplotlib.gridspec as gsp
import matplotlib.patches as mpatch
import matplotlib.pyplot as plt

import fire_an.makeplots.plot_utils as pu
import fire_an.utils.constants_and_units as c

def plotcomp_mass_ion_BH_noBH(fnmass_noBH='', fnion_noBH='',
                              fnmass_BH='', fnion_BH='',
                              outname='', title=None):
    mapn = {'mass_noBH': fnmass_noBH,
            'ion_noBH': fnion_noBH,
            'mass_BH': fnmass_BH,
            'ion_BH': fnion_BH}
    mapkeys = ['mass_noBH', 'mass_BH',
               'ion_noBH', 'ion_BH']
    examplekey_BH = 'mass_BH'
    examplekey_noBH = 'mass_noBH'

    examplekey_ion = 'ion_noBH'
    examplekey_mass = 'mass_noBH'
    maps = {}
    vmins = {}
    vmaxs = {}
    extents = {}
    rvirs = {}
    mvirs = {}
    if title is None:
        filens = {}
        redshifts = {}
    for key in mapn:
        with h5py.File(mapn[key], 'r') as f:
            map = f['map'][:]
            vmin = f['map'].attrs['minfinite']
            vmax = f['map'].attrs['max']

            box_cm = f['Header/inputpars'].attrs['diameter_used_cm']
            cosmopars = {key: val for key, val in \
                        f['Header/inputpars/cosmopars'].attrs.items()}
            #print(cosmopars)
            if 'Rvir_ckpcoverh' in f['Header/inputpars/halodata'].attrs:
                rvir_ckpcoverh = f['Header/inputpars/halodata'].attrs['Rvir_ckpcoverh']
                rvir_pkpc = rvir_ckpcoverh * cosmopars['a'] / cosmopars['h']
            elif 'Rvir_cm' in f['Header/inputpars/halodata'].attrs:
                rvir_cm = f['Header/inputpars/halodata'].attrs['Rvir_cm']
                rvir_pkpc = rvir_cm / (c.cm_per_mpc * 1e-3)
            xax = f['Header/inputpars'].attrs['Axis1']
            yax = f['Header/inputpars'].attrs['Axis2']
            box_pkpc = box_cm / (1e-3 * c.cm_per_mpc)
            extent = (-0.5 * box_pkpc[xax], 0.5 * box_pkpc[xax],
                      -0.5 * box_pkpc[yax], 0.5 * box_pkpc[yax])
            maps[key] = map
            vmins[key] = vmin
            vmaxs[key] = vmax
            extents[key] = extent
            rvirs[key] = rvir_pkpc
            mvirs[key] = f['Header/inputpars/halodata'].attrs['Mvir_g']
            if title is None:
                filen = f['Header/inputpars'].attrs['snapfiles'][0].decode()
                filen = filen.split('/')
                if 'output' in filen:
                    filen = filen[-3]
                else:
                    filen = filen[-2]
                filens[key] = filen
                redshifts[key] = f['Header/inputpars/cosmopars'].attrs['z']
            if key == examplekey_ion:
                pathn = 'Header/inputpars/maptype_args_dict'
                ion_used = f[pathn].attrs['ion'].decode()

    vmin_mass = min([vmins[key] if 'mass' in key else np.inf for key in mapn])
    vmax_mass = max([vmaxs[key] if 'mass' in key else -np.inf for key in mapn])
    cmap_mass = 'viridis'
    mincol_ion = 13.
    vmin_ion = min([vmins[key] if 'ion' in key else np.inf for key in mapn])
    vmax_ion = max([vmaxs[key] if 'ion' in key else -np.inf for key in mapn])
    if vmin_ion < mincol_ion and vmax_ion > mincol_ion:
        cmap_ion = pu.paste_cmaps(['gist_yarg', 'plasma'], 
                                  [vmin_ion, mincol_ion, vmax_ion])
    else:
        cmap_ion = 'plasma'

    ncols = 2
    nrows = 2
    hspace = 0.2
    wspace = 0.35
    panelsize = 2.5
    if title is None:
        tspace = 1.2
    else:
        tspace = 0.5
    caxw = 0.12 * panelsize
    width_ratios = [panelsize] * ncols + [caxw]
    height_ratios = [panelsize] * nrows
    
    width = sum(width_ratios) * (1. + ncols / (ncols + 1.) * wspace)
    height = sum(height_ratios) * (1. + (nrows - 1)/ (nrows) * hspace) + tspace
    tfrac = 1. - tspace / height
    fig = plt.figure(figsize=(width, height))
    grid = gsp.GridSpec(nrows=nrows, ncols=ncols + 1,
                        hspace=hspace, wspace=hspace, 
                        width_ratios=width_ratios,
                        height_ratios=height_ratios,
                        top=tfrac)
    axes = [fig.add_subplot(grid[i // 2, i%ncols]) for i in range(len(mapn))]
    cax_mass = fig.add_subplot(grid[0, 2])
    cax_ion = fig.add_subplot(grid[1, 2])
    fontsize = 12
    
    if title is None:
        if np.allclose(redshifts[examplekey_ion], 
                       [redshifts[key] for key in mapkeys]):
            title = ['redshift: {z:.3f}, ion: {ion}'.format(\
                     z=redshifts[examplekey_ion], ion=ion_used)]
        else:
            print('files have different redshifts: {}'.format(redshifts))
            title = []
        for key, pre in zip([examplekey_noBH, examplekey_BH],
                            ['left (no BH): ', 'right (BH): ']):
            filen = pre + filens[key]
            maxline = 70
            if len(filen) > maxline:
                splitopts = np.where([char == '_' for char in filen])[0]
                origlen = len(filen)
                filen = filen.split('_')
                rest = origlen
                splitlist=[]
                last = 0
                while rest > maxline:
                    if last == 0:
                        start = 0
                    else:
                        start = splitopts[last]
                    cur = last + sum(splitopts[last + 1:] \
                                     <= start + maxline)
                    if cur == last:
                        msg = 'file name splitting for title failed for {}'
                        raise RuntimeError(msg.format('_'.join(filen)))
                    splitlist.append(cur)
                    rest = rest - (splitopts[cur] - splitopts[last])
                    last = cur
                #print(splitopts)
                #print(splitlist)
                #print(origlen)
                for i in splitlist:
                    filen[i] = filen[i] + '\n'
                filen = '_'.join(filen)
            title.append(filen)
        title = '\n'.join(title)

    xlabel = ['X', 'Y', 'Z'][xax] + ' [pkpc]'
    ylabel = ['X', 'Y', 'Z'][yax] + ' [pkpc]'
    fig.suptitle(title, fontsize=fontsize)

    imgs = {}
    for axi, (ax, mapkey) in enumerate(zip(axes, mapkeys)):
        labelx = axi >= len(mapn) - ncols
        labely = axi % ncols == 0
        if labelx:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        if labely:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(axis='both', labelsize=fontsize-2)
        
        if 'ion' in mapkey:
            vmin = vmin_ion
            vmax = vmax_ion
            cmap = cmap_ion
            axtitle = 'ion, '
        elif 'mass' in mapkey:
            vmin = vmin_mass
            vmax = vmax_mass
            cmap = cmap_mass
            axtitle = 'gas, '
        if '_noBH' in mapkey:
            axtitle = axtitle + 'BH off'
        else:
            axtitle = axtitle + 'BH on'
        mvir_logmsun = np.log10(mvirs[mapkey] / c.solar_mass)
        masspart = '$\\mathrm{{M}}_\\mathrm{{vir}}=10^{{{logm_msun:.1f}}}'+\
                   '\\; \\mathrm{{M}}_{{\\odot}}$'
        masspart = masspart.format(logm_msun=mvir_logmsun)
        axtitle = masspart + '\n' + axtitle
        
        imgs[mapkey] = ax.imshow(maps[mapkey].T, origin='lower', 
                                 interpolation='nearest', vmin=vmin,
                                 vmax=vmax, cmap=cmap, extent=extent)

        patches = [mpatch.Circle((0., 0.), rvirs[mapkey])]
        collection = mcol.PatchCollection(patches)
        collection.set(edgecolor=['red'], facecolor='none', linewidth=1.5)
        ax.add_collection(collection)
        if axi == 0:
            ax.text(1.05 * 2**-0.5 * rvirs[mapkey], 
                    1.05 * 2**-0.5 * rvirs[mapkey], 
                    '$R_{\\mathrm{vir}}$',
                    color='red', fontsize=fontsize)
        ax.text(0.02, 0.98, axtitle, fontsize=fontsize - 1, color='red',
                horizontalalignment='left', verticalalignment='top',
                transform=ax.transAxes)

    plt.colorbar(imgs[examplekey_ion], cax=cax_ion, 
                 extend='neither', orientation='vertical')
    cax_ion.set_ylabel('$\\log_{10} \\, \\mathrm{N} \\; [\\mathrm{cm}]^{-2}$',
                        fontsize=fontsize) 
    plt.colorbar(imgs[examplekey_mass], cax=cax_mass, 
                 extend='neither', orientation='vertical') 
    clabel_mass = '$\\log_{10} \\, \\Sigma_{\\mathrm{gas}}' + \
                  ' \\; [\\mathrm{g} \\, \\mathrm{cm}]^{-2}$'
    cax_mass.set_ylabel(clabel_mass, fontsize=fontsize) 
    if outname is not None:
        plt.savefig(outname, bbox_inches='tight')

def plotcomps_mass_ion_BH_noBH(mapset=2):
    outdir = '/Users/nastasha/ciera/projects_lead/fire3_ionabs/maps_BH_noBH/'
    templatetype = 'all'
    if mapset == 2:
        mdir = '/Users/nastasha/ciera/sim_maps/fire/set2_BH_noBH/'
        template_m12i = 'coldens_{qt}_m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_{bh}_gacc31_fa0.5_snap50_shrink-sph-cen_BN98_2rvir_v1.hdf5'
        template_m13h206 = 'coldens_{qt}_m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_{bh}_gacc31_fa0.5_snap50_shrink-sph-cen_BN98_2rvir_v1.hdf5' 
        bh_on = 'sdp1e-4'
        bh_off = 'sdp1e10'
        qt_mass = 'gas-mass'
        qts_ion = ['O6', 'Ne8', 'N5', 'C2', 'Si2', 'Fe2', 'Mg2', 'Mg10']

        templates = {'m12i': template_m12i,
                     'm13h206': template_m13h206}
        labels_tpl = {'m12i': 'm12i_m6e4, MHDCRspec1, fireCR0, crdiffc690, gacc31_fa0.5',
                      'm13h206': 'm13h206_m3e5, MHDCRspec1, fireCR0, crdiffc690, gacc31_fa0.5'}
        title_tpl= 'Gas and {ion} columns, z=0.5, {template}'
        _outname = 'mapcomp_BH_noBH_set2_{template}_snap50_gas_{ion}.pdf'
        tempkeys = list(template.keys()) 
    elif mapset == 3:
        mdir = '/Users/nastasha/ciera/sim_maps/fire/set3_BH_noBH/'
        template_m12i = 'coldens_{{qt}}_m12i_m6e4_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_{{bh}}_gacc31_fa0.5_snap{snap}_shrink-sph-cen_BN98_2rvir_v1.hdf5'
        template_m13h206 = 'coldens_{{qt}}_m13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR0_Oct142021_crdiffc690_{{bh}}_gacc31_fa0.5_snap{snap}_shrink-sph-cen_BN98_2rvir_v1.hdf5' 
        bh_on = 'sdp1e-4'
        bh_off = 'sdp1e10'
        qt_mass = 'gas-mass'
        qts_ion = ['O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        snaps = [45, 49, 51, 60]
        
        templates = {'m12i_snap{}'.format(snap): \
                     template_m12i.format(snap=snap) for snap in snaps}
        templates.update({'m13h206_snap{}'.format(snap): \
                          template_m13h206.format(snap=snap) \
                          for snap in snaps})
        keys_m12i = ['m12i_snap{}'.format(snap) for snap in snaps]
        keys_m13h206 = ['m13h206_snap{}'.format(snap) for snap in snaps]
        labels_tpl = {key: key + ', MHDCRspec1, fireCR0, crdiffc690, gacc31_fa0.5'\
                      for key in keys_m12i}
        labels_tpl.update({key: key + ', MHDCRspec1, fireCR0, crdiffc690, gacc31_fa0.5'\
                           for key in keys_m13h206})              
        title_tpl= 'Gas and {ion} columns, {template}'
        _outname = 'mapcomp_BH_noBH_set3_{template}_gas_{ion}.pdf'
        tempkeys = keys_m12i + keys_m13h206 
    elif mapset == 4:
        mdir = '/Users/nastasha/ciera/sim_maps/fire/set4_BH_noBH/'
        template = 'coldens_{{qt}}_{ic}_m7e3_MHD_fire3_fireBH_Sep182021_hr_crdiffc690_{{bh}}_gacc31_fa0.5_snap{snap}_shrink-sph-cen_BN98_2rvir_v2.hdf5'
        bh_on = 'sdp2e-4'
        bh_off = 'sdp1e10'
        qt_mass = 'gas-mass'
        qts_ion = ['O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        snaps = [500, 258, 186]
        ics = ['m12m', 'm12f']
        
        templates = {'{}_snap{}'.format(ic, snap): \
                     template.format(snap=snap, ic=ic) \
                     for snap in snaps for ic in ics}
        keys = ['{}_snap{}'.format(ic, snap) for snap in snaps for ic in ics]
        labels_tpl = {key: key + ', m7e3, Sep182021_hr_crdiffc690, gacc31_fa0.5'\
                      for key in keys}          
        title_tpl= 'Gas and {ion} columns, {template}'
        _outname = 'mapcomp_BH_noBH_set4_{template}_gas_{ion}.pdf'
        tempkeys = keys
    elif mapset == 5:
        mdir = '/Users/nastasha/ciera/sim_maps/fire/set5_BH_noBH/'
        template_noBH = 'coldens_{{qt}}_{ic}_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5_snap{snap}_shrink-sph-cen_BN98_2rvir_v2.hdf5'
        template_BH = 'coldens_{{qt}}_{ic}_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}_shrink-sph-cen_BN98_2rvir_v2.hdf5'
        bh_on = 'sdp1e-4'
        bh_off = 'sdp1e10'
        qt_mass = 'gas-mass'
        qts_ion = ['O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        snaps = [45, 50, 60]
        ics = ['m13h206', 'm13h007']
        
        templates_BH = {'{}_snap{}'.format(ic, snap): \
                        template_BH.format(snap=snap, ic=ic) \
                        for snap in snaps for ic in ics}
        templates_noBH = {'{}_snap{}'.format(ic, snap): \
                         template_noBH.format(snap=snap, ic=ic) \
                         for snap in snaps for ic in ics}
        keys = ['{}_snap{}'.format(ic, snap) for snap in snaps for ic in ics]   
        title_tpl = None
        templatetype = 'BHnoBH'
        _outname = 'mapcomp_BH_noBH_set5_{template}_gas_{ion}.pdf'
        tempkeys = keys

    elif mapset == 6:
        mdir = '/Users/nastasha/ciera/sim_maps/fire/set6_BH_noBH/'
        template_noBH = 'coldens_{{qt}}_{ic}_m3e5_MHD_fire3_fireBH_Sep182021_crdiffc690_sdp1e10_gacc31_fa0.5_snap{snap}_shrink-sph-cen_BN98_2rvir_v2.hdf5'
        template_BH = 'coldens_{{qt}}_{ic}_m3e5_MHDCRspec1_fire3_fireBH_fireCR1_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000_snap{snap}_shrink-sph-cen_BN98_2rvir_v2.hdf5'
        bh_on = 'sdp1e-4'
        bh_off = 'sdp1e10'
        qt_mass = 'gas-mass'
        qts_ion = ['O6', 'Ne8', 'Mg10', 'N5', 'Mg2']
        snaps = [45, 50]
        ics = ['m13h002']
        
        templates_BH = {'{}_snap{}'.format(ic, snap): \
                        template_BH.format(snap=snap, ic=ic) \
                        for snap in snaps for ic in ics}
        templates_noBH = {'{}_snap{}'.format(ic, snap): \
                         template_noBH.format(snap=snap, ic=ic) \
                         for snap in snaps for ic in ics}
        keys = ['{}_snap{}'.format(ic, snap) for snap in snaps for ic in ics]        
        title_tpl = None
        templatetype = 'BHnoBH'
        _outname = 'mapcomp_BH_noBH_set5_{template}_gas_{ion}.pdf'
        tempkeys = keys

    for templatekey in tempkeys:
        for ion in qts_ion:
            if templatetype == 'all':
                filetemplate = mdir + templates[templatekey]
                fnmass_noBH = filetemplate.format(qt=qt_mass, bh=bh_off)
                fnion_noBH = filetemplate.format(qt=ion, bh=bh_off)
                fnmass_BH = filetemplate.format(qt=qt_mass, bh=bh_on)
                fnion_BH = filetemplate.format(qt=ion, bh=bh_on)
            elif templatetype == 'BHnoBH':
                filetemplate_BH = mdir + templates_BH[templatekey]
                filetemplate_noBH = mdir + templates_noBH[templatekey]
                fnmass_noBH = filetemplate_noBH.format(qt=qt_mass)
                fnion_noBH = filetemplate_noBH.format(qt=ion)
                fnmass_BH = filetemplate_BH.format(qt=qt_mass)
                fnion_BH = filetemplate_BH.format(qt=ion)
            
            outname = outdir + _outname.format(template=templatekey, ion=ion)

            if title_tpl is None:
                title = None
            else:
                title = title_tpl.format(ion=ion, 
                                         template=labels_tpl[templatekey])
                if len(title) > 50:
                    splitopts = np.where([char == ',' for char in title])[0]
                    optchoice = np.argmin(np.abs(splitopts - 0.5 * len(title)))
                    splitind = splitopts[optchoice]
                    title = (title[:splitind + 1]).strip() + '\n' +\
                            (title[splitind + 1:]).strip()
            #print(fnmass_noBH)
            #print(fnion_noBH)
            #print(fnmass_BH)
            #print(fnion_BH)
            #print(outname)
            #print(title)
            plotcomp_mass_ion_BH_noBH(fnmass_noBH=fnmass_noBH, 
                                      fnion_noBH=fnion_noBH,
                                      fnmass_BH=fnmass_BH, 
                                      fnion_BH=fnion_BH,
                                      outname=outname, 
                                      title=title)