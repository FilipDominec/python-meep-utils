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
    def __init__(self, comment="", simtime=50e-12, resolution=4e-6, cells=1, cell_size=50e-6, padding=20e-6, 
            radius=13e-6, wirethick=0):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "SphereArray"    
        self.src_freq, self.src_width = 1000e9, 2000e9    # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (100e9, 2000e9)    # Which frequencies will be saved to disk
        self.pml_thickness = 20e-6

        self.size_x = cell_size 
        self.size_y = cell_size
        self.size_z = cells*cell_size + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cell_size*cells/2)-padding, (cell_size*cells/2)+padding)
        self.cellcenters = np.arange((1-cells)*cell_size/2, cells*cell_size/2, cell_size)

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials
        self.materials = []  

        tio2 = meep_materials.material_TiO2(where=self.where_sphere) 
        if 'LoLoss' in comment:           
            tio2.pol[0]['gamma'] /= 10.   ## optionally edit the first TiO2 optical phonon to have lower damping
        elif 'LossLess' in comment:
            tio2.pol[0]['gamma'] = 0.     ## (or remove damping completely)
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

class SphereCDH_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=50e-12, resolution=4e-6, radius=13e-6, eps2=12., spacing=50e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "Sphere"

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.simtime = simtime      # [s]
        self.src_freq, self.src_width = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = spacing, spacing, spacing

        ## Define materials
        self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        #self.materials = [meep_materials.material_Metal_THz(where = self.where_TiO2)]  
        self.test_materials()

        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        #meep_utils.plot_eps(self.materials, mark_freq=[f_c])


    def where_TiO2(self, r):
        if  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius) and not  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius*.75):
            return self.return_value             # (do not change this line)
        #if  ((in_ycyl(r, cx=0, cz=0, rad=self.radius) and not in_ycyl(r, cx=0, cz=0, rad=self.radius*.75)) and 
                #not (r.z()>0 and in_xslab(r, cx=0, d=self.size_z*.15+self.resolution)) and in_yslab(r, cy= self.size_y*(-.25), d=self.radius*.2)):
            #return self.return_value             # (do not change this line)
        #if  in_yslab(r, cy=self.size_y*(.25), d=self.wtth) and in_zslab(r, cz=0, d=self.wlth):
            #return self.return_value             # (do not change this line)
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
        self.src_freq, self.src_width = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = self.resolution*2, spacing, spacing

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.test_materials()

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
        self.src_freq, self.src_width = 2000e9, 4000e9     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.size_x, self.size_y, self.size_z  = spacing, spacing, zsize

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_TiO2, eps=eps2)]  
        self.materials = [meep_materials.material_Metal_THz(where = self.where_metal)]  
        self.test_materials()

        #f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        #meep_utils.plot_eps(self.materials, mark_freq=[f_c])


    def where_metal(self, r):
        return 0
        if (in_zslab(r, cz=-15e-6, d=self.resolution*2) or in_zslab(r, cz=-15e-6, d=self.resolution*2)) and \
                ((in_yslab(r, cy=0, d=self.radius*2 ) and in_xslab(r, cx=0, d=self.radius*2*self.thin)) or 
                 (in_yslab(r, cy=0, d=self.radius*2*self.thin) and in_xslab(r, cx=0, d=self.radius*2 ))):
            return self.return_value             # (do not change this line)
#}}}
