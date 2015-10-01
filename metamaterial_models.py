#!/usr/bin/env python
#coding:utf8

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

class SphereArray(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=50e-12, resolution=4e-6, cellsize=50e-6, cellnumber=1, padding=20e-6, 
            radius=13e-6, wirethick=0, loss=1, epsilon=-1):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "SphereArray"    
        self.src_freq, self.src_width = 1000e9, 4000e9    # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (100e9, 2000e9)    # Which frequencies will be saved to disk
        self.pml_thickness = 20e-6

        self.size_x = cellsize 
        self.size_y = cellsize
        self.size_z = cellnumber*cellsize + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cellsize*cellnumber/2)-padding, (cellsize*cellnumber/2)+padding)
        self.cellcenters = np.arange((1-cellnumber)*cellsize/2, cellnumber*cellsize/2, cellsize)

        self.register_locals(locals())          ## Remember the parameters
        print self.comment

        ## Define materials (with manual Lorentzian clipping) 
        self.materials = []  

        if epsilon==-1:     ## use titanium dioxide if permittivity not specified...
            tio2 = meep_materials.material_TiO2(where=self.where_sphere) 
            if loss != 1: tio2.pol[0]['gamma'] *= loss   ## optionally modify the first TiO2 optical phonon to have lower damping
        else:           ## ...or define a custom dielectric if permittivity not specified
            tio2 = meep_materials.material_dielectric(where=self.where_sphere, eps=self.epsilon) 

        self.fix_material_stability(tio2, f_c=2e13, verbose=0) ## rm all osc above the first one, to optimize for speed 
        self.materials.append(tio2)

        if wirethick > 0:
            au = meep_materials.material_Au(where=self.where_wire)
            self.fix_material_stability(au, verbose=0)
            self.materials.append(au)

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_sphere(self, r):
        for cellc in self.cellcenters:
            if  in_sphere(r, cx=0, cy=0, cz=cellc, rad=self.radius):
                return self.return_value             # (do not change this line)
        return 0
    def where_wire(self, r):
        for cellc in self.cellcenters:
            if  in_xcyl(r, cy=self.size_y/2, cz=cellc, rad=self.wirethick) or \
                    in_xcyl(r, cy= -self.size_y/2, cz=cellc, rad=self.wirethick):
                return self.return_value             # (do not change this line)
        return 0
#}}}
class RodArray(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=4e-6, cellsize=100e-6, cellnumber=1, padding=20e-6, 
            radius=10e-6, eps2=100):

        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "RodArray"

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.simtime = simtime      # [s]
        self.src_freq, self.src_width = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y  = self.resolution*2, cellsize
        self.size_z = cellnumber*cellsize + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cellsize*cellnumber/2)-padding, (cellsize*cellnumber/2)+padding)
        self.cellcenters = np.arange((1-cellnumber)*cellsize/2, cellnumber*cellsize/2, cellsize)

        ## Define materials
        self.materials = [meep_materials.material_TiO2(where = self.where_TiO2)]  
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  

        for m in self.materials: self.fix_material_stability(m)
        self.test_materials()

    def where_TiO2(self, r):
        #if  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius) and not  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius*.75):
        if  in_xcyl(r, cy=0, cz=0, rad=self.radius):
            return self.return_value             # (do not change this line)
        return 0
#}}}
class Slab(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=2e-6, cellnumber=1, cellsize=100e-6, padding=50e-6,
            fillfraction=0.5, epsilon=2):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "Slab"    
        self.src_freq, self.src_width = 1000e9, 4000e9  # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (10e9, 2000e9)    # Which frequencies will be saved to disk
        self.pml_thickness = 0.1*c/self.src_freq

        self.size_x = resolution*2 
        self.size_y = resolution*2
        self.size_z = cellnumber*cellsize + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cellsize*cellnumber/2)-padding, (cellsize*cellnumber/2)+padding)
        self.cellcenters = np.arange((1-cellnumber)*cellsize/2, cellnumber*cellsize/2, cellsize)

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials
        # note: for optical range, it was good to supply f_c=5e15 to fix_material_stability
        if 'Au' in comment:           
            m = meep_materials.material_Au(where=self.where_metal)
            self.fix_material_stability(m, verbose=0) ## rm all osc above the first one, to optimize for speed 
        elif 'Ag' in comment:           
            m = meep_materials.material_Ag(where=self.where_metal)
            self.fix_material_stability(m, verbose=0) ## rm all osc above the first one, to optimize for speed 
        else:
            m = meep_materials.material_dielectric(where=self.where_metal, loss=.0001)
        self.materials = [m]

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_metal(self, r):
        if in_zslab(r, d=self.cellsize*self.fillfraction, cz=0):
            return self.return_value             # (do not change this line)
        return 0
#}}}
class SRRArray(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=50e-12, resolution=4e-6, cellsize=100e-6, cellnumber=1, padding=20e-6, 
            radius=40e-6, wirethick=6e-6, srrthick=6e-6, splitting=16e-6, splitting2=0e-6, capacitorr=0e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "SRRArray"    
        self.src_freq, self.src_width = 1000e9, 4000e9    # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (100e9, 2000e9)    # Which frequencies will be saved to disk
        self.pml_thickness = 20e-6

        self.size_x = cellsize 
        self.size_y = cellsize
        self.size_z = cellnumber*cellsize + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cellsize*cellnumber/2)-padding, (cellsize*cellnumber/2)+padding)
        self.cellcenters = np.arange((1-cellnumber)*cellsize/2, cellnumber*cellsize/2, cellsize)

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials (with manual Lorentzian clipping) 
        self.materials = []  

        au = meep_materials.material_Au(where=self.where_wire)
        self.fix_material_stability(au, verbose=0)
        self.materials.append(au)

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_wire(self, r):
        for cellc in self.cellcenters:
            ## define the wires
            if  in_xcyl(r, cy=self.size_y/2, cz=cellc, rad=self.wirethick) or \
                    in_xcyl(r, cy= -self.size_y/2, cz=cellc, rad=self.wirethick):
                return self.return_value             # (do not change this line)

        #r = self.RotatedCoordsY(r, angle=np.pi/4)
        #print dir(r)
        r = self.rotatedX(r, np.pi/4)
        r = self.rotatedY(r, np.pi/4)
        #without rot =  66  s
        #with    rot = 165  s
        #with   2rot = 239  s
        #with    rot = 74.7442  s (optimized):
        #with   2rot = 74.7442  s (optimized):

        for cellc in self.cellcenters:
            ## define the split-ring resonator
            if  (((in_ycyl(r, cx=0, cz=cellc, rad=self.radius+self.srrthick/2)             # outer radius
                    and not in_ycyl(r, cx=0, cz=cellc, rad=self.radius-self.srrthick/2)   # subtract inner radius
                    and in_yslab(r, cy=0, d=self.srrthick))                               # delimit to a disc
                    or (in_xcyl(r, cy=0, cz=cellc+self.radius, rad=self.capacitorr) and in_xslab(r, cx=0, d=self.splitting+2*self.wirethick)))
                    and not (r.z()>cellc and in_xslab(r, cx=0, d=self.splitting))           # make the first splitting
                    and not (r.z()<cellc and in_xslab(r, cx=0, d=self.splitting2))):        # make the second splitting, if any
                return self.return_value             # (do not change this line)
        return 0
#}}}
class TMathieu_Grating(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=200e-15, resolution=20e-9, cellnumber=1, padding=50e-6, 
            tdist=50e-6, ldist=100e-6, rcore1=6e-6, rclad1=0, rcore2=6e-6, tshift=0):
        """ I have a red laser (spot size : 2mm of diameter) that goes through 2 grids placed a 50cm (see pictures below) but 100um apart from each other. A photomultiplier is placed behind the grids at 1m. During the experiment the second grid moves transversally and alternatively block the light and let the light reaching the photomultiplier. The grids induce a diffraction pattern, of which we only collect the central bright spot with the photomultiplier (a pinhole is placed in front of it with a 2mm diameter hole). What I would like to do is to simulate the profile of intensity of the light collected at the photomultiplier while the second grid moves. That's a first thing. Secondly, I would like to simulate how the profile changes while some material is deposit on the bars of the first grids and obstruct slowly the light to go through. So what matters to me is to recorded "how much light" of the initial light reach my PM while the second grid moves for various thicknesses of material deposited on the first one. """
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "TMathieu_Grating"    
        self.src_freq, self.src_width = 500e12, 2000e12    # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (380e12, 730e12)    # Which frequencies will be saved to disk
        self.pml_thickness = 500e-9

        self.size_x = resolution*1.8 
        self.size_y = tdist
        self.size_z = ldist + 2*padding + 2*self.pml_thickness
        meep.master_printf("number of voxels: %d", int(self.size_x*self.size_y*self.size_y/resolution**3))
        self.monitor_z1, self.monitor_z2 = (-(ldist/2)-padding, (ldist/2)+padding)
        cellsize = ldist+2*padding

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials (with manual Lorentzian clipping) 
        self.materials = []  

        au = meep_materials.material_Au(where=self.where_wire)
        self.fix_material_stability(au, verbose=0)
        self.materials.append(au)

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_wire(self, r):
        if  in_xcyl(r, cy=0, cz=-self.ldist/2, rad=self.rcore1):                        ## first grid
            return self.return_value             # (do not change this line)

        if  in_xcyl(r, cy=self.tshift, cz=self.ldist/2, rad=self.rcore2) or \
                in_xcyl(r, cy=self.tshift-self.size_y, cz=self.ldist/2, rad=self.rcore2):       ## second grid may be transversally shifted
            return self.return_value             # (do not change this line)

        return 0
#}}}
class HalfSpace(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-15, resolution=10e-9, cellnumber=1, padding=200e-9, cellsize = 200e-9,
            epsilon=33.97, blend=0):
        """ This structure demonstrates that scatter.py can also be used for samples on a substrate with unlimited thickness. The back
        side of the substrate is not simulated, and it is assumed there will be no coherent interferences between its sides.

        To enable the simulations, the monitor planes are enabled to be placed also inside a dielectric. In which case the wave amplitude is 
        adjusted so that the light intensity is maintained. The field amplitudes and phases have physical meaning only when both monitor planes are
        in the same medium, though.

        Besides, the example demonstrates that with the choice of permittivity of ((1+.5**.5)/(1-.5**.5))**2 ~ 33.97 for a steep interface with air, 
        the transmitted and reflected waves have exactly the same energy.
        """
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "HalfSpace"    
        self.src_freq, self.src_width = 500e12, 100e12    # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (10e12, 1000e12)    # Which frequencies will be saved to disk
        self.pml_thickness = 500e-9

        self.size_x = resolution*1.8 
        self.size_y = resolution*1.8
        self.size_z = blend + 2*padding + 2*self.pml_thickness + 6*resolution
        self.monitor_z1, self.monitor_z2 = (-padding, padding)
        self.register_locals(locals())          ## Remember the parameters
        self.mon2eps = epsilon                  ## store what dielectric is the second monitor embedded in

        ## Define materials
        self.materials = []  
        if 'Au' in comment:         self.materials += [meep_materials.material_Au(where=self.where_m)]
        elif 'Ag' in comment:       self.materials += [meep_materials.material_Ag(where=self.where_m)]
        elif 'metal' in comment:    
            self.materials += [meep_materials.material_Au(where=self.where_m)]
            self.materials[-1].pol[1:] = []
            self.materials[-1].pol[0]['gamma'] = 0
        else:                       self.materials += [meep_materials.material_dielectric(where=self.where_m, eps=self.epsilon)]

        for m in self.materials: 
            self.fix_material_stability(m, f_c=3e15) ## rm all osc above the first one, to optimize for speed 

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_m(self, r):
        ## Just half-space
        #if r.z() > 0: return self.return_value

        ## Smooth sine-like transition from air to dielectric: a broadband anti-reflex layer
        if r.z()<-self.blend*.5: return 0
        if r.z()> self.blend*.5 or self.blend==0: return self.return_value
        return self.return_value*(1.+np.sin(r.z()/0.5/self.blend*np.pi/2))/2

        ## Single antireflex layer on substrate
        #if r.z() < 0 and r.z() > -self.padding/2:
            #return self.return_value**.5
        #if r.z() > 0:
            #return self.return_value
        return 0
#}}}

models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray, 'SRRArray':SRRArray, 'TMathieu_Grating':TMathieu_Grating, 'HalfSpace':HalfSpace}

