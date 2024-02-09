import sys

import fire_an.simlists as sl
import fire_an.spectra.genspectra as gs

## Notes on requirements:
# on Quest, 200 sightlines ran in a bit under 6 hours (1 node/1 CPU)
# but only about 300 made it in a 10 hours time limit
#    (2procs/node, 2 CPUs/ proc)
# 16 hours:

def runtest2():
    '''
    first look: spectra around 1 halo
    '''
    simname = 'crheatfix_m12f_r7100'
    snapnum = 277
    obase = f'/test2/tridentray_{simname}_{snapnum}'
    setargs = {'totlength': 4.,
               'gridside': 4.,
               'gridpoints_side': 15, 
               'axis':'z'}
    gs.runsightlines(simname, snapnum, outname_base=obase,
                      settype='grid', **setargs)
    
def runtest4(ind):
    '''
    second look: multiple haloes,
    but only FIRE-2 core, 2 redshifts, 1 axis,
    still using the FG09 UV/X-ray bkg, 
    working on the yt dl fix branch
    '''

    simnames = sl.m12_f2md
    snapshots = [sl.snaps_f2md[0], sl.snaps_f2md[1]]

    snapi = ind % len(snapshots)
    simi = ind // len(snapshots)

    simname = simnames[simi]
    snapnum = snapshots[snapi]
    outbase = f'/test4/tridentray_{simname}_{snapnum}'
    setargs = {'totlength': 4.,
               'gridside': 4.,
               'gridpoints_side': 20, 
               'axis':'z'}
    gs.runsightlines(simname, snapnum, outname_base=outbase,
                      settype='grid', **setargs)
    
def main(ind):
    if ind == 0:
        runtest2()
    elif ind >= 1 and ind < 21:
        runtest4(ind - 1)
    else:
        raise ValueError(f'nothing specified for index {ind}')
    
if __name__ == '__main__':
    ind = int(sys.argv[1])
    main(ind)