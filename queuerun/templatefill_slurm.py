#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
not using the standard .format stuff, since this fills in a bash
scripts and those are just full of ${varname} stuff

Parameters:
-----------
a method specification, e.g. --firemaps_frontera_seq

for --firemaps_frontera_seq: fire_maps.py calls with sequential indices
    --PYFILL_TIME_HHMMSS=HH:MM:SS
    --PYFILL_PARTITION=<flex OR small>
    --PYFILL_JOBNAME=<job name>
    --start=<first index, int>
    --step=<number of tasks per node, int>
    --last=<last index, int>

    generates a set of slurm scripts with tasks divided over them. Each
    calls fire_maps.py with a number of indices, so that all indices 
    from start to last are called, with max. step on each node. No jobs 
    with more than one node are created.

for --frontera_seq_multigen, --stampede2_seq_multigen: 
    fire_maps.py calls with sequential indices
    unlike --firemaps_frontera_seq, the step size determines the number
    of concurrent processes, but the launcher file contains all jobs 
    from start to last. 
    --PYFILL_TIME_HHMMSS=HH:MM:SS
    --PYFILL_PARTITION=<flex OR small>
    --PYFILL_JOBNAME=<job name>
    --start=<first index, int>
    --last=<last index, int>
    --npar=<number of tasks per node, int> (number run in parallel)
'''

import os
import stat
import sys

# where the templates are and the sbatch files go
sdir_frontera = '/work2/08466/tg877653/frontera/slurm/'
sdir_stampede2 = '/work2/08466/tg877653/stampede2/slurm/'

def fillin(templatefilen, outfilen, **kwargs):
    '''
    reads in templatefilen, replaces instances of each kwarg key
    with the corresponding value, writes out the result to outfilen.
    '''
    #print(templatefilen)
    with open(templatefilen, 'r') as f:
        template = f.read()
    out = template
    for key in kwargs:
        val = kwargs[key]
        out = out.replace(key, str(val))
    with open(outfilen, 'w') as f_out:
        f_out.write(out)
    # make executable (for me)
    os.chmod(outfilen, 
             stat.S_IRWXU or stat.S_IRGRP or stat.S_IXGRP or stat.S_IRWXO)

def fillin_firemaps_frontera_seqinds(**kwargs):
    '''
    fill in template_slurm_frontera_autolauncher.sh, checking the key
    values
    '''
    templatefilen = sdir + 'template_slurm_frontera_autolauncher.sh'
    # standard format to document, jobname for batch queue submission
    outfilen = sdir + 'slurm_run_{st}_to_{ls}_{jobname}.sh'
    defaults = {'PYFILL_TIME_HHMMSS': '01:00:00',
                'PYFILL_PARTITION': 'flex'}
    keys_req = ['PYFILL_JOBNAME',
                'PYFILL_NTASKS',
                'PYFILL_FMIND_START']
    keys_opt = list(defaults.keys())
    kwargs_next = defaults.copy()
    for key in kwargs:
        if (key not in keys_req) and (key not in keys_opt):
            msg = 'skipping key {}: not an option for this template'
            print(msg.format(key))
            continue
        kwargs_next.update({key: kwargs[key]})
        if key in keys_req:
            keys_req.remove(key)
    if len(keys_req) > 0:
        print('required keys missing: {}'.format(keys_req))
    
    last = kwargs_next['PYFILL_FMIND_START'] \
           + kwargs_next['PYFILL_NTASKS'] - 1
    outfilen = outfilen.format(st=kwargs_next['PYFILL_FMIND_START'],
                               ls=last,
                               jobname=kwargs_next['PYFILL_JOBNAME'])
    
    fillin(templatefilen, outfilen, **kwargs_next)

def fillset_firemaps_frontera_seqinds(**kwargs):
    keys_req = ['start', 'step', 'last']
    for key in keys_req:
        if key not in kwargs:
            raise ValueError('Argument {} missing'.format(key))
    start = int(kwargs['start'])
    step = int(kwargs['step'])
    last = int(kwargs['last'])
    for key in keys_req:
        del kwargs[key]
    _kwargs_next = kwargs.copy()
    numfiles = (last - start) // step + 1
    if (last - start + 1) % step == 0:
        ntaskss = [step] * numfiles
        starts = [start + step * i for i in range(numfiles)]
    else:
        ntaskss = [step if start + (i + 1) * step <= last
                   else last + 1 - (start + i * step)
                   for i in range(numfiles)]
        starts = [start + sum(ntaskss[:i]) for i in range(numfiles)]
    for ntasks, start in zip(ntaskss, starts):
        kwargs_next = _kwargs_next.copy()
        kwargs_next['PYFILL_NTASKS'] = ntasks
        kwargs_next['PYFILL_FMIND_START'] = start
        fillin_firemaps_frontera_seqinds(**kwargs_next)

def fillin_frontera_seq_multigen(**kwargs):
    '''
    fill in template_slurm_frontera_autolauncher_multigen.sh,
    checking the key values
    '''
    templatefilen = sdir + 'template_slurm_frontera_autolauncher_multigen.sh'
    # standard format to document, jobname for batch queue submission
    outfilen = sdir + 'slurm_run_{st}_to_{ls}_{jobname}.sh'
    defaults = {'PYFILL_TIME_HHMMSS': '01:00:00',
                'PYFILL_PARTITION': 'flex'}
    keys_req = ['PYFILL_JOBNAME',
                'start',
                'last',
                'npar']
    keys_opt = list(defaults.keys())
    kwargs_out = defaults.copy()
    for key in kwargs:
        if (key not in keys_req) and (key not in keys_opt):
            msg = 'skipping key {}: not an option for this template'
            print(msg.format(key))
            continue
        kwargs_out.update({key: kwargs[key]})
        if key in keys_req:
            keys_req.remove(key)
    if len(keys_req) > 0:
        print('required keys missing: {}'.format(keys_req))
    
    kwargs_out['PYFILL_FMIND_START'] = kwargs_out['start']
    kwargs_out['PYFILL_FMIND_LAST'] = kwargs_out['last']
    kwargs_out['PYFILL_NTASKS'] = kwargs_out['npar']

    outfilen = outfilen.format(st=kwargs_out['start'],
                               ls=kwargs_out['last'],
                               jobname=kwargs_out['PYFILL_JOBNAME'])
    del kwargs_out['start']
    del kwargs_out['last']
    del kwargs_out['npar']
    fillin(templatefilen, outfilen, **kwargs_out)

if __name__ == '__main__':
    args = sys.argv[1:]
    methodargs = ['--firemaps_frontera_seq',
                  '--frontera_seq_multigen',
                  '--firemaps_stampede2_seq',
                  '--stampede2_seq_multigen']
    methodset = False
    kwargs = {}
    for arg in args:
        if arg in methodargs:
            if methodset:
                raise ValueError('multiple template options specified')
            methodset = True
            continue
        elif arg.startswith('--'):
            arg = arg[2:]
            key, val = arg.split('=')
        else:
            key, val = arg.split('=')
        kwargs.update({key: val})

    if '--firemaps_frontera_seq' in args:
        sdir = sdir_frontera
        fillset_firemaps_frontera_seqinds(**kwargs)
    elif '--frontera_seq_multigen' in args:
        sdir = sdir_frontera
        fillin_frontera_seq_multigen(**kwargs)
    elif '--stampede2_seq_multigen' in args:
        sdir = sdir_stampede2
        fillin_frontera_seq_multigen(**kwargs)
    elif '--firemaps_stampede2_seq' in args:
        sdir = sdir_stampede2
        fillset_firemaps_frontera_seqinds(**kwargs)
    else: 
        msg = 'No known template specified; options are: {}'
        raise ValueError(msg.format(methodargs))