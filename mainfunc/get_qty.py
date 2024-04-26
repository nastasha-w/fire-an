

import numpy as np
import string
import numbers as num

from fire_an.ionrad.ion_utils import Linetable_PS20, atomw_u_dct, elt_atomw_cgs
import fire_an.mainfunc.coords as coords
import fire_an.mainfunc.haloprop as hp
import fire_an.utils.constants_and_units as c
import fire_an.utils.opts_locs as ol


# tested -> seems to work
# dust on/off, redshifts 1.0, 2.8, Z=0.01, 0.0001
# compared FIRE interpolation to neighboring table values
# tested ions sum to 1: lintable=True -> yes, except for molecules,
#                       dust depletion (high nH, low T, more at higher Z)
#                       lintable=False -> no, in some regions of phase 
#                       space, without good physics reasons
def get_ionfrac(snap, ion, indct=None, table='PS20', simtype='fire',
                ps20depletion=True, lintable=True):
    '''
    Get the fraction of an element in a given ionization state in 
    a given snapshot.

    Parameters:
    -----------
    snap: snapshot reader obect
        exact class depends on the simulation
    ion: str
        ion to get the fraction of. Format e.g. 'o6', 'fe17'
    indct: dict or None
        dictionary containing any of the followign arrays
        'filter': bool, size of arrays returned by the snap reader
                  determines which resolution elements to use.
                  If not repesent, all resolution elements are used.
        If not present, the following values are obtained using snap:
        'logT': temperature in log10 K. 
        'lognH': hydrogen number density in log10 particles / cm**3
        'logZ': metal mass fraction in log10 fraction of total mass (no 
                solar scaling)
    table: {'PS20'}
        Which ionization tables to use.
    simtype: {'fire'}
        What format does the simulation reader class snap have?
    ps20depletion: bool
        Take away a fraction of the ions to account for the fraction of the
        parent element depleted onto dust.
    lintable: bool 
        interpolate the ion balance (and depletion, if applicable) in linear
        space (True), otherwise, it's done in log space (False) 

    Returns:
    --------
        the fraction of the parent element nuclei that are a part of the 
        desired ion
    '''
    if simtype == 'fire':
        readfunc = snap.readarray_emulateEAGLE
        prepath = 'PartType0/'
        redshift = snap.cosmopars.z
    else:
        raise ValueError('invalid simtype option: {}'.format(simtype))
    
    if indct is None:
        indct = {}
    if 'filter' in indct: # gas selection, e.g. a spatial region
        filter = indct['filter']
    else:
        filter = slice(None, None, None)
    if 'logT' in indct: # should be in [log10 K]
        logT = indct['logT']
    else:
        logT = np.log10(readfunc(prepath + 'Temperature')[filter])
        tocgs = snap.toCGS
        if not np.isclose(tocgs, 1.):
            logT += np.log10(tocgs)
    if 'lognH' in indct: # should be in log10 cm**-3
        lognH = indct['lognH']
    else:
        hdens = readfunc(prepath + 'Density')[filter]
        d_tocgs = snap.toCGS
        hmassfrac = readfunc(prepath + 'ElementAbundance/Hydrogen')[filter]
        hmassfrac_tocgs = snap.toCGS
        hdens *= hmassfrac 
        hdens *= d_tocgs * hmassfrac_tocgs / (c.atomw_H * c.u)
        del hmassfrac
        lognH = np.log10(hdens)
        del hdens
    if table in ['PS20']:
        if 'logZ' in indct: # no solar normalization, 
            #just straight mass fraction
            logZ = indct['logZ']
        else:
            logZ = readfunc(prepath + 'Metallicity')[filter]    
            logZ = np.log10(logZ)
            tocgs = snap.toCGS
            if not np.isclose(tocgs, 1.):
                logZ += np.log10(tocgs)
        # Inputting logZ values of -np.inf (zero metallicity, does 
        # happen) leads to NaN ion fractions in interpolation.
        # Since the closest edge of the tabulated values is used anyway
        # it's safe to substute a tiny value like -100.
        if np.any(logZ == -np.inf):
            logZ = logZ.copy()
            logZ[logZ == -np.inf] = -100.
    if table == 'PS20':
        interpdct = {'logT': logT, 'lognH': lognH, 'logZ': logZ}
        iontab = Linetable_PS20(ion, redshift, emission=False, vol=True,
                 ionbalfile=ol.iontab_sylvia_ssh, 
                 emtabfile=ol.emtab_sylvia_ssh, lintable=lintable)
        ionfrac = iontab.find_ionbal(interpdct, log=False)
        if ps20depletion:
            ionfrac *= (1. - iontab.find_depletion(interpdct))
    else:
        raise ValueError('invalid table option: {}'.format(table))
    ## debug map NaN values
    #if np.any(ionfrac) < 0.:
    #    print('some ion fractions < 0 from get_ionfrac')
    #    print('min/max ion fraction values: {}, {}'.format(np.min(ionfrac),
    #                                                       np.max(ionfrac)))
    #if np.any(np.isnan(ionfrac)):
    #    msg = 'Some ion fractions were NaN: {} out of {}'
    #    print(msg.format(np.sum(np.isnan(ionfrac)), len(ionfrac)))
    # if np.any(np.isnan(iontab.iontable_T_Z_nH)):
    #     print('NaN values in the table to be interpolated')
    #     print('Parameters used in Linetable_PS20:')
    #     print('ion: ', ion)
    #    print('redshift: ', redshift)
    #    print('emission: ', False)
    #    print('vol: ', True)
    #    print('ionbalfile: ', ol.iontab_sylvia_ssh)
    #    print('emtabfile: ', ol.emtab_sylvia_ssh)
    #    print('lintable: ', lintable)
    return ionfrac

# untested, including lintable option and consistency with table values
# do a test like test_ionbal_calc before using
# (note that element abundance rescaling and volume multiplication will
# make the direct table comparison harder than with the ion fractions;
# best to do some direct input tests for that)
def get_loglinelum(snap, line, indct=None, table='PS20', simtype='fire',
                   ps20depletion=True, lintable=True, ergs=False,
                   density=False):
    '''
    Get the luminosity (density) of a series of resolution elements.

    Parameters:
    -----------
    snap: snapshot reader obect
        exact class depends on the simulation
    line: str
        line to calculate the luminosity of. Should match the line
        list in the table.
    indct: dict or None
        dictionary containing any of the followign arrays
        'filter': bool, size of arrays returned by the snap reader
                  determines which resolution elements to use.
                  If not repesent, all resolution elements are used.
        If not present, the following values are obtained using snap:
        'logT': temperature in log10 K. 
        'lognH': hydrogen number density in log10 particles / cm**3
        'logZ': metal mass fraction in log10 fraction of total mass (no 
                solar scaling)
        'eltmassf': mass fraction of the line-producing 
                element (no solar scaling)
        'hmassf': log mass fraction of hydrogen (no solar scaling)
        'mass': mass in g
        'density': density in g/cm**3
                
    table: {'PS20'}
        Which ionization tables to use.
    simtype: {'fire'}
        What format does the simulation reader class snap have?
    ps20depletion: bool
        Take away a fraction of the ions to account for the fraction of the
        parent element depleted onto dust.
    lintable: bool 
        interpolate the ion balance (and depletion, if applicable) in linear
        space (True), otherwise, it's done in log space (False) 
    ergs: bool
        output luminosity in erg/s[/cm**3] (True); 
        otherwise, output in photons/s[/cm**3] (False)
    density: bool
        output luminosity density ([erg or photons]/s/cm**3) (True);
        otherwise, output (total) luminosity ([erg or photons]/s) (False)

    Returns:
    --------
        the log luminosity (erg/s or photons/s, depending on ergs value)
        or log luminosity density (erg/s/cm**3 or photons/s/cm**3, depending
        on the ergs value), depending on the density value
        
    '''
    if simtype == 'fire':
        readfunc = snap.readarray_emulateEAGLE
        prepath = 'PartType0/'
        redshift = snap.cosmopars.z
    else:
        raise ValueError('invalid simtype option: {}'.format(simtype))
    
    # read in filter, any arrays already present
    if indct is None:
        indct = {}
    if 'filter' in indct: # gas selection, e.g. a spatial region
        filter = indct['filter']
    else:
        filter = slice(None, None, None)
    if 'logT' in indct: # should be in [log10 K]
        logT = indct['logT']
    else:
        logT = np.log10(readfunc(prepath + 'Temperature')[filter])
        tocgs = snap.toCGS
        if not np.isclose(tocgs, 1.):
            logT += np.log10(tocgs)
    if 'Hmassf' in indct:
        hmassf = indct['Hmassf']
    else:
        hmassf = readfunc(prepath + 'ElementAbundance/Hydrogen')[filter]
        hmassf_tocgs = snap.toCGS
    if 'lognH' in indct: # should be in log10 cm**-3
        lognH = indct['lognH']
    else:
        hdens = readfunc(prepath + 'Density')[filter]
        d_tocgs = snap.toCGS
        hdens *= hmassf
        hdens *= d_tocgs * hmassf_tocgs / (c.atomw_H * c.u)
        del hmassfrac
        lognH = np.log10(hdens)
        del hdens
    if table in ['PS20']:
        if 'logZ' in indct: # no solar normalization, 
            #just straight mass fraction
            logZ = indct['logZ'].copy()
        else:
            logZ = readfunc(prepath + 'Metallicity')[filter]    
            logZ = np.log10(logZ)
            tocgs = snap.toCGS
            if not np.isclose(tocgs, 1.):
                logZ += np.log10(tocgs)
        # interpolation needs finite values, float32(1e-100) == 0
        logZ[logZ == -np.inf] = -100.
    if table == 'PS20':
        interpdct = {'logT': logT, 'lognH': lognH, 'logZ': logZ}
        table = Linetable_PS20(line, redshift, emission=False, vol=True,
                               ionbalfile=ol.iontab_sylvia_ssh, 
                               emtabfile=ol.emtab_sylvia_ssh, 
                               lintable=lintable)
        # log10 erg / s / cm**3 
        luminosity = table.find_logemission(interpdct)
        # luminosity in table = (1 - depletion) * luminosity_if_all_elt_in_gas
        # so divide by depleted fraction to get undepleted emission
        if not ps20depletion:
            luminosity -= \
                np.log10(1. - table.find_depletion(interpdct))
        
        # table values are for solar element ratios at Z
        # rescale to actual element density / hydrogen density
        parentelt = string.capwords(table.element)
        if parentelt == 'Hydrogen':
            del logZ
            del logT
            del lognH
        else:
            linelum_erg_invs_invcm3 -= \
                table.find_assumedabundance(interpdct, log=True)
            del logT
            del lognH
            del logZ

            if 'eltmassf' in indct:
                eltmassf = indct['eltmassf']
            else:
                readpath = prepath + 'ElementAbundance/' + parentelt
                eltmassf = readfunc(readpath)[filter]
            zscale = eltmassf / hmassf
            del eltmassf
            del hmassf
            zscale *= atomw_u_dct['Hydrogen'] / table.elementmass_u
            luminosity += np.log10(zscale)
        if not density:
            # log10 erg / s / cm**3 -> erg / s
            if 'mass' in indct:
                logmass = np.log10(indct['mass'])
                m_toCGS = 1.
            else:
                logmass = np.log10(readfunc(prepath + 'Mass')[filter])
                m_toCGS = snap.toCGS
                        # log10 erg / s / cm**3 -> erg/s
            if 'density' in indct:
                logdens = np.log10(indct['density'])
                d_toCGS = 1.
            else:
                logmass = np.log10(readfunc(prepath + 'Density')[filter])
                d_toCGS = snap.toCGS
            logvol = logmass - logdens
            del logmass
            del logdens
            v_toCGS = np.log10(m_toCGS / d_toCGS)
            if not np.isclose(v_toCGS, 0.):
                logvol += v_toCGS
            luminosity += logvol
            del logvol
        if not ergs:
            # erg -> photons
            wl = table.wavelength_cm
            erg_per_photon = c.planck * c.c / wl
            luminosity -= erg_per_photon         
    else:
        raise ValueError('invalid table option: {}'.format(table))
    return luminosity
    
def get_qty(snap, parttype, maptype, maptype_args, filterdct=None):
    '''
    calculate a quantity to map

    Parameters:
    -----------
    snap: Firesnap object (readin_fire_data.py)
        used to read in what is needed
    parttype: {0, 1, 4, 5}
        particle type
    maptype: {'Mass', 'Volume', 'Metal', 'ion', 'sim-direct', 'coords'}
        what sort of thing are we looking for
    maptype_args: dict or None
        additional arguments for each maptype (dictionary, keys are the
        listed strings)
        for maptype value:
        'Mass': None (ignored)
        'Volume': None (ignored)
            only for gas (parttype 0), as density is required
        'Metal':
            number of nuclei or nucleus density (all ions together)
            'element': str
                element name, e.g. 'oxygen'
            'density': bool
                get the metal number density instead of number of nuclei.
                The default is False.
        'ion':
            number of ions or ion density
            'ion': str
                ion name. format e.g. 'o6', 'fe17'
            'ionfrac-method': {'PS20', 'sim'}. The default is 'PS20'.
                how to calculate the ion fractions
                'PS20': interpolate the Ploeckinger & Schaye (2020) 
                        table
                'sim': read the ion fraction in from the snapshot
            'ps20depletion': bool
                deplete a fraction of the element onto dust and include
                that factor in the ion fraction. Depletion follows the
                Ploeckinger & Schaye (2020) table values.
                The default is False.
                (ignored unless the 'ps20table' calculation is used)
            'lintable': bool
                interpolate the tables in linear space (True) or log 
                space (False). The default is True.
                (ignored unless the 'ps20table' calculation is used)
            'density': bool
                get the ion density instead of number of nuclei.
                The default is False.
        'sim-direct': str
            a quantity stored directly in the simulation, or calculated
            by the simulation snapshot class (e.g., Temperature)
            'field': str
                the name of the field (after 'PartType<#>') to read 
                in from the simulation.
        'coords': (only one of 'pos', 'cen' allowed per item)
            'pos': [0, 1, 2, 'allcart', 'rcen', 'phi', 'azimuth']
                0, 1, 2: position along the axis with this index
                'allcart': for all three of these cartesian axes
                'rcen': distance to the center
                The center is assumed to match the halo center used
                for the map.
            'vel': [0, 1, 2, 'allcart', 'vrad', 'vtot', 'phi',
                    'azimuth']
                 0, 1, 2: velocity along the axis with this index
                'allcart': for all three of these cartesian axes
                'vrad': radial velocity (relative to coordinate
                        center)
                'vtot': total velocity (rms coordinate velocties)
            'multiple': list of dicts
                instead of 'vel' or 'pos', specify a list of 
                dictionaries with one 'vel' or 'pos' key each.
                in this case, the function returns lists of
                values, toCGS, and documenting dict
                (outer list is values or toCGS or dict, inner list
                 matches the different requested coordinates by index)
            'center_cm': length-3 iterable of floats
                Where to center the positions, along the simulation
                x, y, z axes. Required argument.
                Units: physical cm
            'vcen_cmps': length-3 iterable of floats
                Where to center the velocities, along the simulation
                x, y, z axes. Required argument for 'vel' quantities.
                Units: physical cm/s
            'rotmatrix': float array, shape (3, 3)
                matrix to rotate the simulation coordinates. Done after
                centering, before anything else.

    Returns:
    --------
    qty: array
        the desired quantity for all the filtered resolution elements
    toCGS:
        factor to convert the array to CGS units
    todoc:
        dictonary with useful info to store   
        always contains a 'units' entry
        for sim-direct read-in, might just be 'cgs' though
    '''
    basepath = 'PartType{}/'.format(parttype)
    filter = slice(None, None, None)
    todoc = {}
    if filterdct is not None:
        if 'filter' in filterdct:
            filter = filterdct['filter']

    if maptype == 'Mass':
        qty = snap.readarray_emulateEAGLE(basepath + 'Masses')[filter]
        toCGS = snap.toCGS
        todoc['units'] = 'g'
    elif maptype == 'Volume':
        qty = snap.readarray_emulateEAGLE(basepath + 'Masses')[filter]
        toCGS = snap.toCGS
        qty /= snap.readarray_emulateEAGLE(basepath + 'Density')[filter]
        toCGS = toCGS / snap.toCGS
        todoc['units'] = 'cm**3'
        todoc['method'] = 'Masses / Density'
    elif maptype == 'Metal':
        element = maptype_args['element']
        if element == 'total':
            eltpath = basepath + 'Metallicity'
        else:
            eltpath = basepath + 'ElementAbundance/' + string.capwords(element)
        if 'density' in maptype_args:
            output_density = maptype_args['density']
        else:
            output_density = False
        qty = snap.readarray_emulateEAGLE(eltpath)[filter]
        toCGS = snap.toCGS
        if output_density:
            qty *= snap.readarray_emulateEAGLE(basepath + 'Density')[filter]
        else:
            qty *= snap.readarray_emulateEAGLE(basepath + 'Masses')[filter]
        toCGS = toCGS * snap.toCGS
        if element == 'total':
            # number of nuclei would depend on the metal fractions
            # -> that's a pain, just return grams ( / cm^3)
            todoc['units'] = 'g'
        else:
            todoc['units'] = '(# nuclei)'
            toCGS = toCGS / elt_atomw_cgs(element)
        
        if output_density:
            todoc['units'] += ' * cm**-3'
        todoc['density'] = output_density
    elif maptype == 'ion':
        if parttype != 0 :
            msg = 'Can only calculate ion fractions for gas (PartType0),' + \
                   ' not particle type {}'
            raise ValueError(msg.format(parttype))
        ion = maptype_args['ion']
        if 'ionfrac-method' in maptype_args:
            ionfrac_method = maptype_args['ionfrac-method']
        else:
            ionfrac_method = 'PS20'
        if 'density' in maptype_args:
            output_density = maptype_args['density']
        else:
            output_density = False
        simtype = 'fire'
        if ionfrac_method == 'PS20':
            if 'ps20depletion' in maptype_args:
                ps20depletion = maptype_args['ps20depletion']
            else:
                ps20depletion = False
            if 'lintable' in maptype_args:
                lintable = maptype_args['lintable']
            else:
                lintable = True
            # no tables read in here, just an easy way to get parent 
            # element etc.
            dummytab = Linetable_PS20(ion, snap.cosmopars.z, emission=False,
                                      vol=True, lintable=lintable)
            element = dummytab.element
            eltpath = basepath + 'ElementAbundance/' + string.capwords(element)
            qty = snap.readarray_emulateEAGLE(eltpath)[filter]
            toCGS = snap.toCGS
            if output_density:
                dpath = basepath + 'Density'
                qty *= snap.readarray_emulateEAGLE(dpath)[filter]
            else:
                mpath = basepath + 'Masses'
                qty *= snap.readarray_emulateEAGLE(mpath)[filter]
            toCGS =  toCGS * snap.toCGS
            ionfrac = get_ionfrac(snap, ion, indct=filterdct, 
                                  table=ionfrac_method, 
                                  simtype=simtype, ps20depletion=ps20depletion,
                                  lintable=lintable)
            qty *= ionfrac
            toCGS = toCGS / (dummytab.elementmass_u * c.u)
            todoc['table'] = dummytab.ionbalfile
            todoc['tableformat'] = ionfrac_method
            todoc['units'] = '(# ions)'
        if ionfrac_method == 'sim':
            if simtype == 'fire' and ion == 'H1':
                eltpath = basepath + 'ElementAbundance/Hydrogen'
                qty = snap.readarray_emulateEAGLE(eltpath)[filter]
                toCGS = snap.toCGS
                if output_density:
                    dpath = basepath + 'Density'
                    qty *= snap.readarray_emulateEAGLE(dpath)[filter]
                else:
                    mpath = basepath + 'Masses'
                    qty *= snap.readarray_emulateEAGLE(mpath)[filter]
                toCGS = toCGS * snap.toCGS
                hfpath = basepath + 'NeutralHydrogenAbundance'
                qty *= snap.readarray_emulateEAGLE(hfpath)[filter]
                toCGS = toCGS * snap.toCGS
                # just for the element mass
                dummytab = Linetable_PS20(ion, snap.cosmopars.z, 
                                          emission=False,
                                          vol=True, lintable=True)
                toCGS = toCGS / (dummytab.elementmass_u * c.u)
                todoc['info'] = ('neutral H fraction from simulation'
                                 ' NeutralHydrogenAbundance')
                todoc['units'] = '(# ions)'
            else:    
                msg = ('simulation read-in of ion fractions is not available'
                       'for simulation {} and ion {}')
                raise ValueError(msg.format(simtype, ion))
        if output_density:
            todoc['units'] += ' * cm**-3'
        todoc['density'] = output_density
        todoc['ionfrac-method'] = ionfrac_method
    elif maptype == 'sim-direct':
        field = maptype_args['field']
        qty = snap.readarray_emulateEAGLE(basepath + field)[filter]
        toCGS = snap.toCGS
        todoc['units'] = '(cgs {})'.format(field)
        if field == 'Pressure':
            todoc['info'] = 'thermal pressure only'
    elif maptype == 'coords':
        if 'pos' in maptype_args:
            if 'vel' in maptype_args:
                msg = ('specify "vel" XOR "pos" for coordinate quantities. '
                       f'gave maptype_args:\n{maptype_args}')
                raise ValueError(msg)
            indct = [{'pos': maptype_args['pos']}]
            vcen_cmps = None
        elif 'vel' in maptype_args:
            indct = [{'vel': maptype_args['vel']}]
            if 'vcen_cmps' in maptype_args:
                vcen_cmps = maptype_args['vcen_cmps']
            else:
                msg = ('specify a velocity center (vcen_cmps) for velocities.'
                       f' gave maptype_args:\n{maptype_args}')
                raise ValueError(msg)
        elif 'multiple' in maptype_args:
            indct = []
            vcen_set = False
            for cdct in maptype_args['multiple']:
                if len(cdct) != 1:
                    msg = ('dictionaries in the "multiple" list should only'
                           ' contain a single "pos" or "vel" entry. This '
                           'argument was given as:'
                           f'\n{maptype_args["multiple"]}')
                    raise ValueError(msg)
                if 'vel' in cdct:
                    if 'vcen_cmps' in maptype_args:
                        vcen_cmps = maptype_args['vcen_cmps']
                        vcen_set = True
                    else:
                        msg = ('specify a velocity center (vcen_cmps)'
                               ' for velocities. gave maptype_args:'
                               f'\n{maptype_args}')
                        raise ValueError(msg)
                indct.append(cdct)
            if not vcen_set:
                vcen_cmps = None
        else: 
            msg = ('No coordinate type(s) specified ("pos", "vel", '
                   'or "multiple") in maptype_args. gave maptype_args:\n'
                   f'{maptype_args}')
            raise ValueError(msg)
        if 'center_cm' in maptype_args:
            center_cm = maptype_args['center_cm']
        else:
            msg = ('specify a center (center_cm) for coordinate quantities. '
                   f'gave maptype_args:\n{maptype_args}')
            raise ValueError(msg)
        if 'rotmatrix' in maptype_args:
            rotmatrix = maptype_args['rotmatrix']
        else:
            rotmatrix = None
        coordobj = coords.CoordinateWranger(snap, center_cm, 
                                            rotmatrix=rotmatrix,
                                            parttype=parttype, periodic=False, 
                                            vcen_cmps=vcen_cmps,
                                            filterdct=filterdct)
        out = coordobj.calccoords(indct)
        if len(indct) == 1:
            qty, toCGS, todoc = out[0]
        else:
            qty = [out[i][0] for i in range(len(indct))]
            toCGS = [out[i][1] for i in range(len(indct))]
            todoc = [out[i][2] for i in range(len(indct))]
    else:
        raise ValueError('Invalid maptype: {}'.format(maptype))
    return qty, toCGS, todoc

def process_typeargs_coords(simpath, snapnum, typeargs,
                            paxis=None):
    '''
    allows standard values or gethalodata_shrinkinssphere/get_vcom
    arguments to be used instead of specifiying vcen_cmps, rcen_cm
    explicitly.

    Parameters
    ----------
    typeargs: dict
        may match any maptype_args for maptype 'coords' in get_qty
        additional convenience options (by key):
        'rcen_cm': not included or None
            not included or None: use the BN98 overdensity definition,
            shrinking-sphere method with all particle types (except 2)
            and standard parameter values for centering.
            Otherwise, should match the get_qty options.
        'vcen_cmps': not included, None, or dict
            centering/virial radius always match the defaults under 
            'rcen_cm', and the virial radius is always the 'BN98' one.
            The defaults for 'radius_rvir' and 'parttypes' are 1. and
            'all', respectively. If a dict argument is given, these
            defaults may be overridden by including the 'radius_rvir'
            and/or 'parttypes' keywords, with options matching those
            arguments for haloprop.get_vcom 
            Otherwise, should match the get_qty options.
        'pos': 'los'
            get the line-of-sight position along whatever the 
            projection axis is. (requires the paxis argument to be set
            to 0, 1, or 2)
        'vel': 'los', 'doplos'
            get the line-of-sight velocity (or doppler velocity) 
            along whatever the 
            projection axis is. (requires the paxis argument to be set
            to 0, 1, or 2)
    paxis: {0, 1, 2, None}
        the projection axis. Required if a 'pos' or 'vel' coordinate
        is specified as 'los', otherwise ignored.
    '''
    needsv = 'vel' in typeargs or \
             ('multiple' in typeargs 
              and np.any(['vel' in dct for dct in typeargs['multiple']]))
    outdoc = {}
    typeargs_out = typeargs.copy()
    if ('center_cm' not in typeargs) or \
       ('center_cm' in typeargs and typeargs['center_cm'] is None):
        rhdat, rdoc = hp.gethalodata_shrinkingsphere(simpath, snapnum, 
                                                     meandef='BN98')
        outdoc.update({'coords_' + key: rdoc[key] for key in rdoc})
        rcen_cm = np.array([rhdat['Xc_cm'], rhdat['Yc_cm'], rhdat['Zc_cm']])
        typeargs_out.update({'center_cm': rcen_cm})
        outdoc.update({'coords_rcen_cm_in': 'default'})
        outdoc.update({'coords_center': 'shrinksph'})
    if needsv:
        if 'vcen_cmps' not in typeargs:
            typeargs['vcen_cmps'] = dict()
        elif typeargs['vcen_cmps'] is None:
            typeargs['vcen_cmps'] = dict()
        if (isinstance(typeargs['vcen_cmps'], dict)):
            tvdct = typeargs['vcen_cmps'].copy()
            if 'radius_rvir' in tvdct:
                radius_rvir = tvdct['radius_rvir']
            else:
                radius_rvir = 1.
            outdoc.update({'vcen_radius_rvir': radius_rvir})
            if 'parttypes' in tvdct:
                parttypes = tvdct['parttypes']
            else:
                parttypes = 'all'
            outdoc.update({'vcen_parttypes': parttypes})
            vdat, vdoc = hp.get_vcom(simpath, snapnum, 
                                     radius_rvir, meandef_rvir='BN98',
                                     parttypes=parttypes)
            outdoc.update({'vcen_' + key: vdoc[key] for key in vdoc})
            vcen_cmps = np.array([vdat['VXcom_cmps'],
                                  vdat['VYcom_cmps'],
                                  vdat['VZcom_cmps']])
            typeargs_out.update({'vcen_cmps': vcen_cmps})
        else:
            vcin = typeargs['vcen_cmps']
            if not (hasattr(vcin, '__len__') 
                    and len(vcin) == 3
                    and np.all([isinstance(vcin[i], num.Number) 
                                for i in range(3)])):
                raise ValueError('The "vcen_cmps" argument should be a'
                                 ' length 3 iterable of floats, '
                                 'a dictionary, or None')
    if 'vel' in typeargs and typeargs['vel'] == 'los':
        outdoc.update({'coords_vel_in': 'los'})
        typeargs_out['vel'] = paxis
    elif 'vel' in typeargs and typeargs['vel'] == 'doplos':
        outdoc.update({'coords_veldop_in': 'los'})
        typeargs_out['vel'] = f'dop{paxis}'
    elif 'pos' in typeargs and typeargs['pos'] == 'los':
        outdoc.update({'coords_pos_in': 'los'})
        typeargs_out['pos'] = paxis
    elif 'multiple' in typeargs:
        cspec = typeargs['multiple'].copy()
        for di, dct in enumerate(cspec):
            if 'vel' in dct and dct['vel'] == 'los':
                outdoc.update({'coords_vel_in': 'los'})
                cspec[di].update({'vel': paxis})
            elif 'vel' in dct and dct['vel'] == 'doplos':
                outdoc.update({'coords_veldop_in': 'los'})
                cspec[di].update({'vel': f'dop{paxis}'})
            elif 'pos' in dct and dct['pos'] == 'los':
                outdoc.update({'coords_pos_in': 'los'})
                cspec[di].update({'pos': paxis})
        typeargs_out['multiple'] = cspec
    return typeargs_out, outdoc