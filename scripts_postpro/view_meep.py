#!/usr/bin/env python2.7
#-*- coding: utf-8 -*-
import os,sys
import numpy as np

try: from enthought.mayavi import mlab
except: from mayavi import mlab
mlab.options.offscreen = True ## XXX

def symmetric_colors(lut_manager):
    """ Given a lookup table, sets its colormap to symmertric red-white-blue. Useful to view signed data."""
    lut_manager.use_default_range = False
    max_data_range = max(abs(lut_manager.data_range))
    lut_manager.data_range = np.array([-max_data_range, max_data_range])
    #lut_manager.lut_mode = lut_mode
def color_black(lut_manager):
    lut_manager.use_default_range = False
    lut_manager.data_range = np.array([0,0])

def RunMayavi(directory):
    figure = mlab.figure(bgcolor=(.5,.5,.5))
    #figure.scene.disable_render = True
    #figure.scene.magnification = 5
    #figure.scene.off_screen_rendering = True

    ## Try to load the structure as wireframe
    if os.path.exists(os.path.join(directory,'structure.vtk')):
        data_source = mlab.pipeline.open(os.path.join(directory,'structure.vtk'))
        iso_surface = mlab.pipeline.iso_surface(data_source, color=(.8,.8,.8), contours=2, opacity=.3)
        iso_surface.actor.property.representation = 'wireframe'
        #iso_surface.actor.property.representation = 'surface'
        iso_surface.actor.property.line_width = .5
        #iso_surface.contour.minimum_contour = 2.
        iso_surface.name = 'Permittivity 3-D contours'

        surface = mlab.pipeline.surface(data_source, opacity=.3, colormap='Greys', vmin=1.0, vmax=10.)

        ## For black outline in 2-D view
        #surface.actor.property.specular_color = (0.0, 0.0, 0.0)
        #surface.actor.property.diffuse_color = (0.0, 0.0, 0.0)
        #surface.actor.property.ambient_color = (0.0, 0.0, 0.0)
        #surface.actor.property.color = (0.0, 0.0, 0.0)
        #surface.actor.property.representation = 'wireframe'

        surface.name = 'Permittivity surface'
        surface.visible = False

    mlab.orientation_axes(opacity=.5)
    mlab.title(os.path.basename(directory), size=.2, height=.01)






    ## Process each file
    filenames = [n for n in os.listdir(directory) if (n[-4:]=='.vtk' and not ('structure' in n) and (not 'slice' in n))]
    filenames.sort()
    for (count, filename) in enumerate(filenames):
        ## Load the files
        data_source = mlab.pipeline.open(os.path.join(directory, filename))

        ## Plot electric field as the x-component scalar value (blue-white-red)
        if 'Evec' in filename:
            evs = mlab.pipeline.extract_vector_components(data_source)
            cutplane = mlab.pipeline.scalar_cut_plane(evs, 
                    plane_orientation='x_axes', colormap='RdBu')
            symmetric_colors(evs.children[0].scalar_lut_manager)
            cutplane.actor.mapper.interpolate_scalars_before_mapping = True
            cutplane.implicit_plane.widget.enabled = False
            cutplane.name = "Electric field x-component"

            ## Warp scalar 
            cutplane.enable_warp_scalar = True
            cutplane.warp_scalar.filter.normal = np.array([ 1.,  0.,  0.])
            cutplane.warp_scalar.filter.scale_factor = 100.0


        ## Plot magnetic field as black vectors
        if 'Hvec' in filename:
            #cutplane = mlab.pipeline.vector_cut_plane(data_source, 
                    #plane_orientation='x_axes', mask_points=1, scale_factor=2, mode="cone", color=(0,0,0))
            cutplane = mlab.pipeline.vector_cut_plane(data_source, 
                    plane_orientation='x_axes', mask_points=1, scale_factor=2, mode="cone")
            cutplane.name = "Magnetic field"
            cutplane.implicit_plane.widget.enabled = False

        if count>1: cutplane.visible = False        ## hide all snapshots except the first

    ## Set up the camera
    figure.scene.camera.clipping_range = [1, 1000]
    figure.scene.camera.focal_point = [10, 10, 30]

    # Y-Z view 
    figure.scene.camera.position = [140, 10, 27]
    figure.scene.camera.focal_point = [10, 10, 27]
    figure.scene.camera.view_up = [0, -1, 0]

    # X-Z view 
    #figure.scene.camera.position = [10, 140, 27]
    #figure.scene.camera.view_up = [1, 0, 0]
    #figure.scene.camera.compute_view_plane_normal()

    figure.scene.parallel_projection = True



    figure.scene.disable_render = False
    figure.scene.render()
    figure.scene.save(u'/home/filip/snapshot4.png')

    mlab.show()

if __name__ == '__main__':		# avoid execution if loaded as a module only
    cwd = os.getcwd()
    if len(sys.argv)>1: 
        directory = sys.argv[1]
    elif os.path.exists(os.path.join(cwd, 'last_simulation_name.txt')):
        #print "Loading from", os.path.join(cwd, 'last_simulation_name.txt')
        directory = os.path.join(cwd, open(os.path.join(cwd, 'last_simulation_name.txt'),'r').read())
    else:
        print "Defaulting to cwd"
        directory = cwd
    RunMayavi(directory)
