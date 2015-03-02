## Introduction
MEEP is an implementation of the FDTD algorithm, which simulates the propagation of electromagnetic
waves. During my PhD studies, I used it intensively for simulations of metamaterials, photonic crystals and
many other problems.

Unlike many commercial packages, the simulation has to be *programmed*; MEEP is only a library accessible from C/C++, Scheme or Python. I chose to write my simulations using Python, as it has powerful modules  
for scientific computing. Starting with MEEP was not easy, though, mostly due to lack of resources online. Over few 
years I wrote multiple utility functions that in my opinion greatly facilitate the simulations, and I publish 
them under open license, with the hope they might be useful to the scientific community.

I will be happy if these scripts help you with your thesis, paper or just any project. In such a case, you can 
made a reference to my website or send me a message; perhaps I can even help you with some useful tips.

Filip Dominec, filip.dominec@gmail.com,
2012 - 2015


## File overview
#### General modules and other files
 * `meep_utils.py`       - the main module with routines useful for python-meep simulations
 * `meep_materials.py`   - module containing realistic definition of materials used 
 * `README.md`		     - this file
 * `LICENSE`			 - GPLv2
 * `plot_scan_as_contours.py` - if multiple simulations are run as a parametric scan, this allows to present the results in a single contour plot
 * `effparam.py`		 - retrieves the metamaterial effective parameters from the complex reflection and transmission (e.g. from `scatter.py`)
 * `plot_cdh.py`,`plot_cdh_new.py` - plots data for current-driven homogenization, TODO fix 

#### Simulation scripts
 * `scatter.py`			 - defines a metamaterial cell containing a dielectric sphere, and optionally metallic wires parallel to electric field
 * `cdh.py` - TODO
 * `spdc.py` - TODO

#### Examples using the simulation scripts
Usually, everything you need to run an example is to change to its directory, and launch `./batch.sh`. In a multiprocessing environment, it is recommended to launch it like `mpirun -np 4 ./batch.sh`. 

 * [x]  `example_metamaterial_s_parameters/` - computes effective parameters of a metamaterial (using `scatter.py` and `effparam.py`); shows how the negative index of refraction is achieved by adding wires, and how it retains/changes when more metamaterial cells are computed (which however can suffer from wrong branch switch)
 * [x]  `example_frequency_domain_solver/` - runs `scatter.py` multiple times in frequency-domain, and then compares the results to the classical Fourier-transformed time-domain simulation
 * [x]  `example_cylindrical_cavity_modes/` - defines a metallic cavity
 * [ ]  `example_surface_plasmons/`
	TODO add support for metal/diel substrate, and try to show the sym-asym interference
 * [ ]  `example_aperture_near-field_microscope/` 
    TODO
 * [ ]  `example_dielectric_bars_width_scan/`
	TODO
 * [ ]  `example_dielectric_slab_oblique_incidence/`
	TODO , c.f. transfer-matrix
 * [ ]  `example_refraction_on_MM_wedge_2D/` - defines a wedge of a 2-D rod array (studied earlier both as a photonic crystal and a metamaterial), and by the means of spatial Fourier transform, analyzes how a beam is refracted depending on its frequency. Compares the result with the s-parameter retrieval method.
	TODO implement seamless 2-D support
 * [ ]  `example_nonlinear_Kerr_focusing/` - demonstrates a source with custom spatial shape, which launches a focused Gaussian beam. Different amplitudes are scanned to show how the nonlinear medium changes the beam and eventually allows filament formation.
	TODO implement nonlinearity, test out
 * [ ]  `example_SPDC/`
	TODO

## Related resources
 * Official website of MEEP: http://ab-initio.mit.edu/wiki/index.php/Meep     
Information on the FDTD algorithm and simulations in general, documentation of the MEEP functions. Examples 
are mostly in Scheme. 
   * Since 2014, the MEEP source code is hosted at Github: https://github.com/stevengj/meep
   * Your questions may (or may not) be answered at the mailing list: meep-discuss@ab-initio.mit.edu,      
http://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/

 * Website of the python-meep interface: https://launchpad.net/python-meep     
Provides examples of how the python-meep functions can be used in scripts.

 * I also write own website on simulations: http://f.dominec.eu/meep/index.html     
My experience with installation requirements and procedure, simulation performance, realistic definition of 
materials, data postprocessing etc. is elaborated there.

 * License: GPLv2, http://www.gnu.org/licenses/gpl-2.0.html


## TODO
- [x] scatter.py, cdh and others should output sim_param in the header (moreover CDH has weird header!!)
- [x] move Kx, Ky out of the model parameters
- [x] put the models into separate module
- [ ] sync harminv from its module with meep_utils, and remove from the latter
- [ ] why I do not see interference of sym/asym plasmons in the example? wrong metal model?

- [ ] from scipy.misc import imsave; imsave('../docs/static/tutorial-epsilon.png', -N.rot90(epsilon)) ? 
- [ ] Use average_field_function instead of my own averaging!
- [ ] use synchronize_fields() instead of shifting H(t) ? - benchmark
- [ ] test averaging on SRR
- [ ] test the Fresnel inversion algorithm on dispersive dielectric slabs
- [ ] fix the stupid SWIG bug: http://sourceforge.net/donate/?user_id=246059#recognition
- [x] resonant modes extraction via HarmInv, done in a branched file 
- [ ] optimize the structure using D.E (http://inspyred.github.com) or CMA-ES
- [ ] mode separation on the user-defined ports
- [x] add examples (tests / case study?):
   * waveguide-splitter
   * metamaterial parameters of dielectric rods (CASE STUDY)
   * metamaterial parameters of dielectric sphere in wire mesh (CASE STUDY)
   * a split-ring resonator and current-driven homogenisation
   * surface-plasmons
   * surface-plasmons on thin-metal (CASE STUDY) 
   * thin-gold-film-transmission
   * plasmonic resonance in gold nanoparticles
   * resistive-metal strips
   * extraordinary transmission
   * Kerr nonlinearity and self-focusing
   * scattering SNOM microscope (CASE STUDY)
   * oblique-wave fabry-pérot resonances, comparison with analytic solution
   * resonances in cylinder cavity, application of harminv and comparison with analytic
   * modeling spontaneous parametric down-conversion
