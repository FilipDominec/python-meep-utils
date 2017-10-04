# Example: angle-frequency scan
This simulation example runs a series of time-domain simulations of the same periodic structure (e.g. an interface of two materials, multiple layers, a diffraction grating, or a layer of a metamaterial). The simulations differ by the transverse wavenumber that is computed: When this wavenumber is zero, we obtain a simulation with perpendicular incidence like most other examples. 

When the transverse wavenumber *K* is nonzero, the source amplitude is spatially modulated as a complex exponential, that is exp(i *K x*) = sin(*K x*) + i cos(*K x*). This corresponds to a plane wave being emitted under an oblique angle, which is determined simply by the Pythagorean theorem as arcsin[*K* / (ω/*cn*)], where *ω* and *n* are the frequency of the wave and the refractive index of the medium surrounding the source.

Obviously, in a time-domain simulation where a broadband pulse is emitted, the angle is *frequency dependent*. At the high-frequency portion of the spectrum, the wave propagates close to the *z*-axis (almost perpendicular to the source plane), whereas at moderate frequencies it deviates much more. At very low frequencies, where *K* > (ω/*cn*), the wave is not radiated at all! This is all true, and it is not a problem if one is interested in a full 2-D scan over a band of frequencies and over a wide variety of angles. One only needs to position each data point at the correct angle, and interpolate the contours.

# Usage
Change, e.g. to the ```01_diel_eps09_blend0u``` directory and run the `./batch.sh` script. After about 45 simulations (each with a different ```Kx``` parameter), one big 2-D map of angle-frequency scan is plotted.


## TODO: Optical pulse reflecting from a graduated interface
In the definition of the structure, the python-meep-utils scripts enable one to *blend* the materials, i.e. specify the density between 0 and 100 % of the material.

This can be used to build gradient lenses, custom-shaped perfectly absorbing objects, model arbitrary doping profiles in semiconductors etc.

A particularly interesting example is a continously blended boundary with a plasma-like medium, illuminated by a plane wane. 
When the plasma is lossy (e.g. gold or silver in the optical range), the wave is continuously damped and only a little part reflects. However, if the plasma is defined by a lossless Drude model, the wave penetrates some depth where the group velocity becomes zero, then the energy is accumulated at some point and eventually radiated back.

The time-coordinate record can be visualized using Paraview:
![Time-coordinate plot of the optical pulse reflected from a gradient plasma wall](./plasma_wall.png)

Note that the field amplitude is about 2.5x stronger at the peak, than was that of the illuminating wave. 

Note also how the phase velocity increases as the relative permittivity drops from 1 towards zero, which manifests itself as the time-coordinate plot being curved.

*More details can be found in meep_metamaterials.py in the HalfSpace model.*

