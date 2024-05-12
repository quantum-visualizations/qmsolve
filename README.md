# QMsolve: A module for solving and visualizing the Schrödinger equation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11181978.svg)](https://doi.org/10.5281/zenodo.11181978)
[![PyPi Python Versions](https://img.shields.io/badge/python-3-1abc9c.svg)](https://pypi.python.org/pypi/qmsolve/)
[![PyPI Version](https://img.shields.io/pypi/v/qmsolve.svg)](https://pypi.python.org/pypi/qmsolve)
[![Conda Version](https://img.shields.io/conda/v/conda-forge/qmsolve)](https://anaconda.org/conda-forge/qmsolve)
[![Chat](https://img.shields.io/static/v1?logo=discord&label=chat&message=on%20discord&color=7289da)](https://discord.gg/xVgFfe7jQ9)

![animation](https://github.com/quantum-visualizations/qmsolve/blob/main/images/3D_two_gaussian_wells.gif)

QMsolve seeks to provide a solid and easy to use solver, capable of solving the Schrödinger equation for one and two particles, 
and creating descriptive and stunning visualizations of its solutions both in 1D, 2D, and 3D.


## Installation

```
pip install qmsolve
```

3D plotting requires to have installed [Mayavi](https://docs.enthought.com/mayavi/mayavi/installation.html). To install Mayavi directly along with QMsolve, you can type:

```
pip install qmsolve[with_mayavi]
```

## Usage

The way this simulator works is by discretizing the Hamiltonian with an arbitrary potential, 
that you can specify as a function of the particle observables. This is achieved with the `Hamiltonian` constructor.

For example:
```python
#define the interaction potential
def harmonic_oscillator(particle):
    k = 100. * eV / Å**2
    return 0.5 * k * particle.x**2

#define the Hamiltonian
H = Hamiltonian(particles = SingleParticle(), 
                potential = harmonic_oscillator, 
                spatial_ndim = 1, N = 512, extent = 20*Å)
```

Then, the method `Hamiltonian.solve` can be used to efficiently diagonalize the Hamiltonian and output the energies and the eigenstates of the system:

```python
eigenstates = H.solve(max_states = 30)
print(eigenstates.energies) # The printed energies are in eV.
```

Finally, the eigenstates can be plotted with the use of the `visualization` class.

The `visualization.superpositions` method features the possibility of interactively visualizing a superposition of the computed eigenstates and studying the time dependence of the resulting wavefunction. 

For efficiently solving the time dependent Schrödinger equation, `TimeSimulation` class must be used. It takes as an argument the Hamiltonian you have previously defined, and the method you desire to use.

For a quick start, take a look to the examples found in the [examples subdirectory](https://github.com/quantum-visualizations/qmsolve/tree/main/examples).



## Eigenstate Solver Examples

To perform the example simulations, just run from the corresponding Python example scripts on the command prompt.


|<img src="https://github.com/quantum-visualizations/qmsolve/blob/main/images/1D_harmonic_oscillator_eigenstates.gif" width="100%">|<img src="https://github.com/quantum-visualizations/qmsolve/blob/main/images/1D_harmonic_oscillator.gif" width="82%">
|:--------------------:|:--------------------:|
`python 1D_harmonic_oscillator.py` | `python 1D_harmonic_oscillator_superpositions.py` |
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/eigenstate%20solver%20examples/1D_harmonic_oscillator.py)| [Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/eigenstate%20solver%20examples/1D_harmonic_oscillator_superpositions.py)|

This is the simplest example and one of the most well-studied Hamiltonians. The first script, `1D_harmonic_oscillator.py` serves as an example of how to use the interface of the Eigenstate solver. This script returns the energies and a visualization of the eigenstates of the harmonic oscillator with an interactive slider. The second script, `1D_harmonic_oscillator_superpositions.py` uses exactly the same Hamiltonian, but this time returns a simulation of a superposition of the computed eigenstates, whose coefficients can be interactively modified using the circular widgets presented below the animation.
 

---
|<img src="https://github.com/quantum-visualizations/qmsolve/blob/main/images/1D_interactive_fermions.gif" width="100%">|<img src="https://github.com/quantum-visualizations/qmsolve/blob/main/images/1D_non_interactive_fermions.gif" width="95%">|
:--------------------:|:--------------------:|
`python 1D_interactive_fermions_HO.py` | `python 1D_non_interactive_fermions_HO.py` |
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/eigenstate%20solver%20examples/1D_interactive_fermions_HO.py)| [Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/eigenstate%20solver%20examples/1D_non_interactive_fermions_HO.py)|

These two examples show two fermions confined in 1D harmonic oscillator. The top-left plot represents the configuration space of the system. 

Because we have two 1D particles, we need a two-dimensional space to represent it. Notice that because the particles are fermions,  the wavefunction satisfies the antisymmetric condition: `𝜓(x1, x2) = - 𝜓(x2, x1)` 

An interesting characteristic that can be observed in these examples is how in the non interactive case the energy levels are equally spaced and degenerated, while in the interactive case the degeneracy is broken.
As a starting point we suggest you to modify the confinement and the interaction potential to see what happens!

The time dependent version of this example can be found in `1D_interactive_fermions_HO_superpositions.py`

---
![animation](https://github.com/quantum-visualizations/qmsolve/blob/main/images/3D_four_gaussian_wells.gif)|
:-----------------------:|
`python 3D_four_gaussian_wells.py`|
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/eigenstate%20solver%20examples/3D_four_gaussian_wells.py)|

This example serves to illustrate how to use the 3D solver. The results for other number of Gaussian wells are presented in [this video](https://www.youtube.com/watch?v=eCk8aIIEZSg). Gaussian wells are preferred as didactic examples because the absence of a singularities in its potential makes the computations easier. For realistic atomic examples, you can take a look at `3D_hydrogen_atom.py` and `3D_dihydrogen_cation.py` which use Coulomb potential.

Furthermore, in the hydrogen atom example you can turn on an electric field to visualize the [Stark effect](https://en.wikipedia.org/wiki/Stark_effect), or a magnetic field in  `3D_hydrogen_atom_magnetic_field.py` to visualize the [Zeeman effect](https://en.wikipedia.org/wiki/Zeeman_effect).

---

The default numerical method used for diagonalizing the Hamiltonian is [ARPACK implementation of Implicitly Restarted Lanczos Method](https://www.caam.rice.edu/software/ARPACK/), which is called with scipy via [eigsh method](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html). For 3D examples, [LOBPCG](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html#rbbbc6164e7a5-4) is preferred since its convergence is faster for high-dimensional matrices.

3D examples are also considerably faster when using a GPU. GPU acceleration requires having [CuPy](https://docs.cupy.dev/en/stable/install.html) and [CUDA](https://developer.nvidia.com/cuda-downloads) installed in your computer. 

To use GPU acceleration in your 3D simulations, add the argument `method ='lobpcg-cupy'` in the Hamiltonian `solve` method. For example:

```python
eigenstates = H.solve( max_states = 50, method ='lobpcg-cupy')
```

## Time Dependent Examples

These examples use the `TimeSimulation` class. This class takes as arguments the Hamiltonian you have previously defined, and the method you desire to use. Currently, there are two methods implemented: [Split-Step Fourier](https://en.wikipedia.org/wiki/Split-step_method) and [Cayley-Crank-Nicolson](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method). To specify the method you want to use, use the arguments `method ='split-step'` and `method ='crank-nicolson'` respectively.

Once the `TimeSimulation` is set up, you need to use `TimeSimulation.run` method to perform the simulation. It takes the following arguments: 

 - `initial_wavefunction`: State of the system at `t=0`, specified as a function of the particle observables.
 - `dt`: size of simulation time step.
 - `total_time`: amount of time to evolve the wavefunction.
 - `store_steps`: number of time steps to save in the array denoted by `TimeSimulation.Ψ` to later visualize or analyze. 
 `total_time`/`store_steps` must be larger or equal to `dt`.

By default in the examples, `initial_wavefunction` initializes a Gaussian wave-packet with a spatial standard deviation equal to `σ` and initial momentum `p_x0` and `p_y0`.


![animation](https://github.com/quantum-visualizations/qmsolve/blob/main/images/2D_double-slit.gif)|
:-----------------------:|
`python 2D_double_slit.py`|
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/time%20dependent%20solver%20examples/2D_double_slit.py)|

This is a famous example, which was used to demonstrate the wave-like behavior of matter. In this script, the user can vary the slits separation, width, and depth to study their effect on the diffracted wavefunction.


---

![animation](https://github.com/quantum-visualizations/qmsolve/blob/main/images/2D_cyclotron_orbit_magneticfield.gif)|
:-----------------------:|
`python 2D_cyclotron_orbit_magneticfield.py`|
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/time%20dependent%20solver%20examples/2D_cyclotron_orbit_magneticfield.py)|


This script shows an example where Crank Nicolson method is required. It presents a charged particle under a strong and uniform magnetic field, being confined in a quantum mechanical cyclotron orbit. The period and the radius of the orbit are compared with the classical values. Unlike other confinement potentials like the harmonic oscillator, the initial wavepacket is greatly dispersed over time.

---


![animation](https://github.com/quantum-visualizations/qmsolve/blob/main/images/2D_quantum_resonant_tunneling.gif)|
:-----------------------:|
`python 2D_quantum_resonant_tunneling.py`|
[Link to the example](https://github.com/quantum-visualizations/qmsolve/blob/main/examples/time%20dependent%20solver%20examples/2D_quantum_resonant_tunneling.py)|

Finally, we present a demonstration of quantum resonant tunneling. In this script, two wavepackets are incident on a double potential well, representing two independent electrons. Despite having less energy than the lower, the upper electron has a higher chance of passing through the barriers. This is because its mean energy has been tuned to excite the quasi-ground state of the double potential well. 

For electrons with energy corresponding approximately to the resonant energy level of the double potential well, the transmittance is close to unity. Furthemore, we can get the energy transmittance spectrum by taking the Fourier transform of the simulated wavepackets at a specific output position, yielding the following plot:

![N|Solid](https://github.com/quantum-visualizations/qmsolve/blob/main/images/quantum_resonant_tunneling.png)
---

Generally, the Split Step method is preferred over Crank Nicolson both because of the computational cost of the time step and its associated error. Split Step has a time step error of cubic order `O(Δt^3)` while Crank Nicolson time step error is quadratic `O(Δt^2)`. Thus Crank Nicolson requires more time steps to achieve the Split Step accuracy. However, Split Step can only be used for potentials of the form `V(x,y,z)` and Crank Nicolson use is necessary when the potential is also dependent of the particles momentum, like in the cyclotron orbit example.

Both methods are also implemented to be GPU-accelerated with cupy. Specifically, the speed boost of the cupy split-step is tested to be one and two orders of magnitudes faster than the CPU implementation, depending of the GPU and the grid size used. To use GPU acceleration in your simulations, use the arguments `method ='split-step-cupy'` and `method ='crank-nicolson-cupy'` in the `TimeSimulation` constructor.

The interface uses [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) for input. In the file [constants.py](https://github.com/quantum-visualizations/qmsolve/blob/main/qmsolve/util/constants.py) there is a list of common conversion factors from other units, that can be imported and used to build your potential.

## Main developers

- [Rafael de la Fuente](https://github.com/rafael-fuente)
- [marl0ny](https://github.com/marl0ny)

## Contributors

- [Hudson Smith](https://github.com/dhudsmith) (shift-invert trick for eigsh solver, lobpcg prototype, support for the project and multiple suggestions)
- [Andrew Knyazev](https://github.com/lobpcg) (lobpcg algorithm and AMG preconditioner)

