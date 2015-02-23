#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """


import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class SphereCDH_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, radius=30e-6, eps2=12., spacing=100e-6, wlth=20e-6, wtth=20e-6, Kx=0, Ky=0, Kz=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "srrstripe"

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = spacing, spacing, spacing

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.materials = [meep_materials.material_Metal_THz(where = self.where_TiO2)]  
        self.TestMaterials()

        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        #meep_utils.plot_eps(self.materials, mark_freq=[f_c])


    def where_TiO2(self, r):
        #if  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius) and not  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius*.75):
        if  ((in_ycyl(r, cx=0, cz=0, rad=self.radius) and not in_ycyl(r, cx=0, cz=0, rad=self.radius*.75)) and 
                not (r.z()>0 and in_xslab(r, cx=0, d=self.size_z*.15+self.resolution)) and in_yslab(r, cy= self.size_y*(-.25), d=self.radius*.2)):
            return self.return_value             # (do not change this line)
        if  in_yslab(r, cy=self.size_y*(.25), d=self.wtth) and in_zslab(r, cz=0, d=self.wlth):
            return self.return_value             # (do not change this line)
        return 0
#}}}
class RodCDH_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, radius=10e-6, spacing=100e-6, eps2=100, Kx=0, Ky=0, Kz=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "srrstripe"

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = self.resolution*2, spacing, spacing

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.TestMaterials()

        #f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        #meep_utils.plot_eps(self.materials, mark_freq=[f_c])


    def where_TiO2(self, r):
        #if  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius) and not  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius*.75):
        if  in_xcyl(r, cy=0, cz=0, rad=self.radius):
            return self.return_value             # (do not change this line)
        return 0
#}}}
class FishnetCDH_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=6e-6, radius=30e-6, spacing=100e-6, zsize=50e-6, thin=.3, Kx=0, Ky=0, Kz=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "fishnet"

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = spacing, spacing, zsize

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.materials = [meep_materials.material_Metal_THz(where = self.where_metal)]  
        self.TestMaterials()

        #f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        #meep_utils.plot_eps(self.materials, mark_freq=[f_c])


    def where_metal(self, r):
        return 0
        if (in_zslab(r, cz=-15e-6, d=self.resolution*2) or in_zslab(r, cz=-15e-6, d=self.resolution*2)) and \
                ((in_yslab(r, cy=0, d=self.radius*2 ) and in_xslab(r, cx=0, d=self.radius*2*self.thin)) or 
                 (in_yslab(r, cy=0, d=self.radius*2*self.thin) and in_xslab(r, cx=0, d=self.radius*2 ))):
            return self.return_value             # (do not change this line)
#}}}

# Model selection
#model = SphereCDH_model(**model_param)
#model = RodCDH_model(**model_param)
model = FishnetCDH_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes="None")

## Create the fields object
f = meep.fields(s)
# Define the Bloch-periodic boundaries (any transversal component of k-vector is allowed)
f.use_bloch(meep.X, -model.Kx/(2*np.pi))
f.use_bloch(meep.Y, -model.Ky/(2*np.pi))
f.use_bloch(meep.Z, -model.Kz/(2*np.pi))

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    src_time_type = meep.band_src_time(model.srcFreq/c, model.srcWidth/c, model.simtime*c/10)
    #src_time_type = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2),
        meep.vec( model.size_x/2,  model.size_y/2, model.size_z/2))
#f.add_volume_source(meep.Ex, src_time_type, srcvolume)

class AmplitudeFactor(meep.Callback): 
    def __init__(self, Kx=0, Ky=0, Kz=0): 
        meep.Callback.__init__(self)
        (self.Kx, self.Ky, self.Kz) = Kx, Ky, Kz
    def complex_vec(self, vec):   ## Note: the 'vec' coordinates are _relative_ to the source center
        ## Current-driven homogenisation source forces the K-vector in whole unit cell
        return np.exp(-1j*(self.Kx*vec.x() + self.Ky*vec.y() + self.Kz*vec.z())) 
af = AmplitudeFactor(Kx=model.Kx, Ky=model.Ky, Kz=model.Kz)
meep.set_AMPL_Callback(af.__disown__())
f.add_volume_source(meep.Ex, src_time_type, srcvolume, meep.AMPL)


## Define monitors planes
#monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':model.Kx, 'Ky':model.Ky}
#monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
#monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)

## Define volume monitors for CDH
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'size_z':model.size_z, 'Kx':model.Kx, 'Ky':model.Ky, 'Kz':model.Kz}
monitor1_Ex = meep_utils.AmplitudeMonitorVolume(comp=meep.Ex, **monitor_options)

slice_makers = [] # [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
#slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_x=0, min_timestep=.05e-12, outputhdf=True, outputpng=True)]
#slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_t=3.6e-12, at_z=model.radius+model.resolution*2, min_timestep=.05e-12, outputhdf=True, outputpng=True)]
#slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_x=0, at_y=0, min_timestep=.03e-12, outputhdf=True, outputpng=True)]
#slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_t=2.5e-12)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step()
    dt = (f.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt)
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):                               # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        for monitor in (monitor1_Ex,): monitor.record(field=f)
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for monitor in (monitor1_Ex,): monitor.record(field=f)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
if meep.my_rank() == 0:
    if  not os.path.exists('cdh'): os.mkdir('cdh')
    with open('cdh/output_Kz=%.3e.dat' % model.Kz, 'w') as outfile: 
        outfile.write("#Parameters Parameters\n")
        outfile.write("#param Kz,%.3e\n" % model.Kz)
        outfile.write("#param simulation_orig_name,%s\n" % model.simulation_name)
        outfile.write("#x-column Frequency [Hz]\n#Column Ex real\n#Column Ex imag\n")
        t,E = monitor1_Ex.get_waveforms()
        np.savetxt(outfile, zip(t, E.real, E.imag), fmt="%.8e")

    #freq, s11, s12 = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            #frequency_domain=sim_param['frequency_domain'], frequency=sim_param['frequency'], 
            #maxf=model.srcFreq+model.srcWidth, pad_zeros=1.0, Kx=model.Kx, Ky=model.Ky)
    #meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
    with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 
    #import effparam        # process effective parameters for metamaterials

meep.all_wait()         # Wait until all file operations are finished
