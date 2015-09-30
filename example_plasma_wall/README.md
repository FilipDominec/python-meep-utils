# Optical pulse reflecting from 
In the definition of the structure, the python-meep-utils scripts enable one to *blend* the materials, i.e. specify the density between 0 and 100 % of the material.

This can be used to build gradient lenses, custom-shaped perfectly absorbing objects, model arbitrary doping profiles in semiconductors etc.

A particularly interesting example is a gradient boundary with a plasma-like medium, illuminated by a plane wane. 
When the plasma is lossy (e.g. gold or silver in the optical range), the wave is continuously damped and nothing reflects. However, if the plasma is defined by a lossless Drude model, the energy penetrates up to some depth where the group velocity becomes zero, the energy is accumulated at some point and eventually radiated back.

The time-coordinate record can be visualized using Paraview:
![Time-coordinate plot of the optical pulse reflected from a gradient plasma wall](./plasma_wall.png)

Note that the field amplitude is about 2.5x stronger at the peak, than was that of the illuminating wave. 

Note also how the phase velocity increases as the relative permittivity drops from 1 towards zero, which manifests itself as the time-coordinate plot being curved.

*More details can be found in meep_metamaterials.py in the HalfSpace model.*

## Usage
Run the `./batch.sh` script
