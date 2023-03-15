import mainfunc.makemap as mm 

# hard to do a true test, but check that projected masses and centering
# sort of make sense
def tryout_massmap(opt=1, center='AHFsmooth'):
    outdir = 'ls'
    _outfilen = 'mass_pt{pt}_{sc}_snap{sn}_ahf-cen_2rvir_v1.hdf5'
    if opt == 1:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 600
    elif opt == 2:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 399
    elif opt == 3:
        parttypes = [0, 1, 4]
        dirpath = '/projects/b1026/snapshots/metal_diffusion/m12i_res7100/'
        simcode = 'metal-diffusion-m12i-res7100'
        snapnum = 196

    for pt in parttypes:
        outfilen = outdir + _outfilen.format(pt=pt, sc=simcode, 
                                             sn=snapnum)
        mm.massmap(dirpath, snapnum, radius_rvir=2., particle_type=pt,
                   pixsize_pkpc=3., axis='z', outfilen=outfilen,
                   center=center)
        
def tryout_wholezoom(index):
    outdir = '/projects/b1026/nastasha/tests/start_fire/map_tests/'

    if index == 0:
        dirpath = '/projects/b1026/snapshots/fire3/m13h206_m3e5/' + \
               'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000/' 
        simname = 'm13h206_m3e5__' + \
                  'm13h206_m3e5_MHDCRspec1_fire3_fireBH_fireCR1' + \
               '_Oct252021_crdiffc1_sdp1e-4_gacc31_fa0.5_fcr1e-3_vw3000'                     
        snapnum = 27  
        outfilen_template = 'mass_pt{pt}_{sc}_snap{sn}_axis-{ax}_' + \
                            'wholezoom_v1.hdf5'
        _temp = outdir + outfilen_template 
        outfilens = {'outfilen_gas': _temp.format(pt=0, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_DM': _temp.format(pt=1, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_stars': _temp.format(pt=4, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),
                     'outfilen_BH': _temp.format(pt=5, sc=simname, 
                                                 sn=snapnum, ax='{ax}'),                            
                    }

    mm.massmap_wholezoom(dirpath, snapnum, pixsize_pkpc=3.,
                      **outfilens)