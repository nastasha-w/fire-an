import sys

import genspectra as gs

def runtest2():
    simname = 'crheatfix_m12f_r7100'
    snapnum = 277
    obase = f'/test2/tridentray_{simname}_{snapnum}'
    setargs = {'totlength': 4.,
               'gridside': 4.,
               'gridpoints_side': 15, 
               'axis':'z'}
    gs.runsightlines(simname, snapnum, outname_base=obase,
                      settype='grid', **setargs)
    
def main(ind):
    if ind == 0:
        runtest2()
    else:
        raise ValueError(f'nothing specified for index {ind}')
    
if __name__ == '__main__':
    ind = int(sys.argv[1])
    main(ind)