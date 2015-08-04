## About this example
MEEP allows to run a FDTD simulation in the time domain, but also has a built-in BiCGStab solver that optimizes the fields shape to match the eigenfunction problem:

Lψ = exp(-iωt) ψ, 

where L is the Maxwell operator [1] ψ is the (vector) field function sought for and ω is the user-defined angular frequency.

Here, the script `scatter.py` is ran multiple times in frequency-domain, and its results are computed to the classical Fourier-transformed time-domain simulation. A big effort was put into the Python scripts to make the switch between FDTD and FDFD simulations as easy as possible, that is, by adding the parameter `frequency=...`.

### References
[1] John D Joannopoulos et al. *Photonic crystals: molding the flow of light*. Princeton university press, 2011

## Usage
Run `batch.sh`. Total computation time is about 20 minutes in single process. 

See its short code for inspiration how to use the frequency-domain simulations.


## Expected results
The results are in the form of one big plot with the comparison, accompanied by the `png/` directory containing the Ex field shapes from the frequency-domain.

![The result of the batch.sh script](example_frequency_domain_solver/SphereArray_simtime=1.000e-10_wirethick=1.000e-05.png)

![The Ex field amplitude at the frequency of 1 THz](./example_frequency_domain_solver/png/At1.000e+12Hz_at_x0.000e+00_tinf.png)
