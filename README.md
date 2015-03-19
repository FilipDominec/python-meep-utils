## Introduction
MEEP is an implementation of the FDTD algorithm, which simulates the propagation of electromagnetic waves. During my PhD studies, I used it intensively for simulations of metamaterials, photonic crystals and many other problems.

Unlike many commercial packages, the simulation has to be *programmed*; MEEP is only a library accessible from C/C++, Scheme or Python. I chose the python-meep, as it enables seamless integration with powerful modules for scientific computing. Starting with python-meep was not easy, though, mostly due to lack of resources online. Over few years I wrote many different simulations, and noticed that they have much code in common. I therefore collected them in a module *meep_utils.py*, that in my opinion greatly facilitates the simulations. To help others starting with python-meep, I also addded several example scripts that are listed below.

I will be happy if these scripts help you with your thesis, paper or just any project. In such a case, you can made a reference to my website or send me a message; perhaps I can even help you with some useful tips.  

Filip Dominec, filip.dominec@gmail.com,
2012 - 2015


## File overview
#### General modules and other files
 * `meep_utils.py`       - the main module with routines useful for python-meep simulations
 * `meep_materials.py`   - module containing realistic definition of materials used 
 * `README.md`		 - this file
 * `LICENSE`		 - General Public License
 * `metamaterial_models.py` - different metamaterial models (that can be shared by other scripts)
 * `plot_scan_as_contours.py` - if multiple simulations are run as a parametric scan, this allows to present all results in a single contour plot
 * `harminv_wrapper.py` - allows to simply use filter diagonalisation method from Python 

#### Examples using the simulation scripts
Usually, everything you need to run an example is to change to its directory, and launch `./batch.sh`. In a multiprocessing environment, it is recommended to launch it like `mpirun -np 4 ./batch.sh`. 

 * [x]  `example_metamaterial_s_parameters/`, `scatter.py`, `effparam.py` - retrieves the effective parameters of a metamaterial using the s-parameters method; shows how the negative index of refraction is achieved by adding wires to a dielectric sphere. Finally plots how the parameters retain/change when more metamaterial cells are stacked.
 * [x]  `example_frequency_domain_solver/` - runs `scatter.py` multiple times in frequency-domain, and then compares the results to the classical Fourier-transformed time-domain simulation
 * [x]  `example_current_driven_homogenisation/`, `cdh.py`, `plot_cdh.py` - computes and plots data for current-driven homogenization; compares them with those obtained from s-parameters
 * [x]  `example_cylindrical_cavity_modes/` - defines a metallic cylindrical cavity, excites the field by a short pulsed source, and analyzes the ringdown to search for all modes. A comparison of Fourier transform, filter-diagonalisation method and the textbook analytic solution is plotted. 
   * [ ]	TODO Subsequently, the shapes of the lowest modes are found by frequency-domain solver.
   * [ ]	fix `../harminv_wrapper.py`, rm harminv from `meep_utils.py`
 * [x]  `example_surface_plasmons/`, `plasmons.py` - an aperture in a thin metal sheet couples incident light to surface plasmons. If the film is surrounded by two media with similar index of refraction, circular interference pattern can be observed between the symmetric and antisymmetric plasmon modes.  A different (hyperbolic) interference pattern can be obtained when the plasmons are coupled by two holes.
   * [ ]	TODO add support for metal/diel substrate, 
   * [ ] and try to show the sym-asym interference
 * [ ]  `example_dielectric_bars_width_scan/`
   * [ ]	TODO
 * [ ]  `example_dielectric_slab_oblique_incidence/`
   * [ ]	TODO , c.f. transfer-matrix
 * [ ]  `example_refraction_on_MM_wedge_2D/` - defines a wedge of a 2-D rod array (studied earlier both as a photonic crystal and a metamaterial), and by the means of spatial Fourier transform, analyzes how a beam is refracted depending on its frequency. Compares the result with the s-parameter retrieval method.
   * [ ]	TODO implement seamless 2-D support
 * [ ]  `example_nonlinear_Kerr_focusing/` - demonstrates a source with custom spatial shape, which launches a focused Gaussian beam. Different amplitudes are scanned to show how the nonlinear medium changes the beam and eventually allows filament formation.
   * [ ]	TODO implement nonlinearity, test out
 * [ ]  `example_SPDC/`, `spdc.py` - TODO
   * [ ]	TODO
 * [ ]  `example_aperture_near-field_microscope/` 
   * [ ] TODO

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
- [ ] effparam.py does not cope with "plot_freq_max=None" anymore? -- fix
- [ ] why I do not see interference of sym/asym plasmons in the example? wrong metal model?
- [ ] plot_contour to read any column from direct sim output / effparam
- [ ] fix 'ValueError: width and height must each be below 32768'
- [ ] stability of metals - try to increase 'gamma' until it goes unstable; map the parameters

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
   * metamaterial parameters of dielectric sphere in wire mesh 
   * a split-ring resonator and current-driven homogenisation
   * surface-plasmons
   * surface-plasmons on thin-metal 
   * thin-gold-film-transmission
   * plasmonic resonance in gold nanoparticles
   * resistive-metal strips
   * extraordinary transmission
   * Kerr nonlinearity and self-focusing
   * scattering SNOM microscope 
   * oblique-wave fabry-pérot resonances, comparison with analytic solution
   * resonances in cylinder cavity, application of harminv and comparison with analytic
   * modeling spontaneous parametric down-conversion
