#!/usr/bin/env python
#-*- coding: utf-8 -*-
import meep_utils, meep_materials
import numpy as np
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab
"""
"""

class Wedge_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=3e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            radius=10e-6, yspacing=100e-6, zspacing=100e-6, monzd=200e-6, epsilon=100):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "Wedge"    ## 
        self.register_locals(locals())          ## Remember the parameters

        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'DielLossless' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.0)]
        elif 'DielLoLoss' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.005)]
        else:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0)]


        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 50e-6
        self.size_x = resolution/1.8
        self.size_y = 2200e-6
        self.size_z = 2000e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.srcWidth = 3000e9     
        self.srcFreq = 4e9 + self.srcWidth/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 2000e9) 

        #meep_utils.plot_eps(self.materials, freq_range=(1e10, 1e14), plot_conductivity=True)
        self.TestMaterials() ## catches (most) errors in structure syntax, before they crash the callback

    def where_wire(self, r):
        y,z = r.y(), r.z()
        if z < self.size_z*(.3) and z-y/2>self.size_z*(-.2): 
            #return self.return_value
            yy,zz = np.arccos(np.cos(y/self.yspacing*np.pi*2))*self.yspacing/(np.pi*2), np.arccos(np.cos(z/self.zspacing*np.pi*2))*self.zspacing/(np.pi*2)
            if (yy**2+zz**2)**.5 < self.radius:
                return self.return_value
        return 0
#}}}
