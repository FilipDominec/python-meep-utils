#!/usr/bin/env python
#coding:utf8    ## důležité: určuje kódování souboru; nutno zvolit cp1250 nebo utf8 dle editoru

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
    def __init__(self, comment="", simtime=15e-15, resolution=4e-9, cellnumber=1, cellsize=50e-9, padding=20e-9):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "Slab"    
        self.src_freq, self.src_width = 1000e12, 4000e12  # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (100e12, 2000e12)    # Which frequencies will be saved to disk
        self.pml_thickness = 20e-9

        self.size_x = resolution*2 
        self.size_y = resolution*2
        self.size_z = cellnumber*cellsize + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cellsize*cellnumber/2)-padding, (cellsize*cellnumber/2)+padding)
        self.cellcenters = np.arange((1-cellnumber)*cellsize/2, cellnumber*cellsize/2, cellsize)

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials
        self.materials = []  
        if 'Au' in comment:           
             self.materials += [meep_materials.material_Au(where=self.where_metal)]
        elif 'Ag' in comment:           
             self.materials += [meep_materials.material_Ag(where=self.where_metal)]
        else:
             self.materials += [meep_materials.material_Ag(where=self.where_metal)]

        for m in self.materials: 
            self.fix_material_stability(m, f_c=5e15, verbose=0) ## rm all osc above the first one, to optimize for speed 

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_metal(self, r):
        if in_zslab(r, d=self.cellsize, cz=0):
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

models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray, 'SRRArray':SRRArray}

