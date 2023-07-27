import os

import fire_an.simlists as sl

simnames_all = sl.m12_hr_all2 + sl.m12_sr_all2 + sl.m12_f2md + \
               sl.m13_hr_all2 + sl.m13_sr_all2

pars_toread = ['Omega0', 'OmegaLambda', 'Omega0', 'HubbleParam',
               'SeedBlackHoleMass',
               'BAL_f_accretion', 
               'BAL_v_outflow', 
               'SeedBlackHolePerUnitMass',
               'BlackHoleRadiativeEfficiency',
               'BH_CosmicRay_Injection_Efficiency',
               'CosmicRay_SNeFraction',
               'UnitVelocity_In_CGS', 'UnitMass_In_CGS']

def get_paramfile(simname):
    simpath = sl.dirpath_from_simname(simname)
    opts_pardir = ['', 'output/']
    opts_parfile = ['params.txt-usedvalues',
                    'parameters-usedvalues',  
                    #'gizmo_parameters.txt-usedvalues',
                    'params.txt',
                    #'gizmo_parameters.txt',
                    ]
    parameterfile = None
    for subdir in opts_pardir:
        for opt in opts_parfile:
            if os.path.isfile(simpath + subdir + opt):
                parameterfile = simpath + subdir + opt
                break 
    if parameterfile is None:
        print(f'Did not find a parameter file for {simname} in {simpath}')
    return parameterfile

def getpars(simname):
    parfilen = get_paramfile(simname)
    if parfilen is None:
       return None
    with open(parfilen, 'r') as f:
        contents = f.readlines()
    lines = contents.split()
    out = {par: None for par in pars_toread}
    for par in pars_toread:
        for line in lines:
            if line.startswith(par):
                val = line.split()[-1]
                try:
                    val = int(val)
                except ValueError:
                    val = float(val)
                out[par] = val
    return out

def saveinfo(outname):
    heads = ['simname'] + pars_toread
    hed = '\t'.join(heads) + '\n'
    linetemp = '\t'.join([f'{{{col}}}' for col in heads]) + '\n'
    with open(outname, 'w') as fo:
        fo.write(hed)
        for simname in simnames_all:
            pars = getpars(simname)
            if pars is None:
                continue
            pars['simname'] = simname
            line = linetemp.format(**pars)
            fo.write(line)




    