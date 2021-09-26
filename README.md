# QMsolve: A module for solving and visualizing the Schrödinger equation
[![PyPi Python Versions](https://img.shields.io/badge/python-3-1abc9c.svg)](https://pypi.python.org/pypi/qmsolve/)
[![PyPI Version](https://img.shields.io/pypi/v/qmsolve.svg)](https://pypi.python.org/pypi/qmsolve)
[![Chat](https://img.shields.io/badge/chat-on%20discord-7289da.svg)](https://discord.gg/xVgFfe7jQ9)

![animation](/images/3D_two_gaussian_wells.gif)

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

## How it works

The way this simulator works is by discretizing the Hamiltonian with an arbitrary potential, 
specified as a function of the particle observables. This is achieved with the `Hamiltonian` constructor.

Then, the `Hamiltonian.solve` method efficiently diagonalizes the Hamiltonian and outputs the energies and the eigenstates of the system.
Finally, the eigenstates can be plotted with the use of the `visualization` class.

The `visualization.superpositions` method features the possibility of interactively visualizing a superposition of the computed eigenstates and studying the time dependence of the resulting wavefunction. 

For efficiently solving the time dependent Schrödinger equation, `TimeSimulation` class must be used. It takes as an argument the Hamiltonian you previuosly built, and the method you desire to use.

For a quick start, take a look to the examples found in the [examples subdirectory](https://github.com/quantum-visualizations/qmsolve/tree/main/examples).

## Eigenstate Solver Examples

To perform the simulations, just run from the corresponding Python scripts on the command prompt.

```
python 1D_harmonic_oscillator.py
```
This is the simplest example and one of the most well-studied Hamiltonians. The script uploaded serves as an example of how to use the interface of the Eigenstate solver. The first part of the script returns the energies and the visualization of the eigenstates of the harmonic oscillator. The second part returns a simulation of a superposition of the computed eigenstates, whose coefficients can be interactively modified using the circular widgets presented below the animation.
 
![animation](/images/1D_harmonic_oscillator.gif)

```
python 1D_interactive_fermions_HO.py
```

These two examples shows two fermions confined in 1D harmonic oscillator. The top-left plot represents the configuration space of the system. Because we have two 1D particles, we need a two-dimensional space to represent it. Notice that because the particles are fermions,  the wavefunction satisfies the antisymmetric condition: `𝜓(x1, x2) = - 𝜓(x2, x1)`
An interesting characteristic you can notice in these examples is how in the non interactive case the energy levels are equally spaced and degenerated, however in the interactive case the degeneracy is broken.
As a starting point we suggest you to modify the confinement and the interaction potential to see what happens!
The time dependent version of this example can be found in `1D_interactive_fermions_HO_superpositions.py`


![animation](/images/1D_interactive_fermions.gif)

```
python 1D_non_interactive_fermions_HO.py
```

![animation](/images/1D_non_interactive_fermions.gif)

```
python 3D_four_gaussian_wells.py
```
This example serves to illustrate how to use the 3D solver. Gaussian wells are prefered as didactic examples because the absence of a singularity in their potential makes the computations easier. For a realistic examples of atoms, you can take a look at `3D_hydrogen_atom.py` and `3D_dihydrogen_cation.py` which use Coulomb potential. 
Furthermore, in the hydrogen atom example you can turn on an electric field to visualize the [Stark effect](https://en.wikipedia.org/wiki/Stark_effect), or a magnetic field in  `3D_hydrogen_atom_magnetic_field.py` to visualize the [Zeeman effect](https://en.wikipedia.org/wiki/Zeeman_effect)

![animation](/images/3D_four_gaussian_wells.gif)

The default numerical method used for diagonalizing the Hamiltonian is [ARPACK implementation of Implicitly Restarted Lanczos Method](https://www.caam.rice.edu/software/ARPACK/), which is called with scipy via [eigsh method](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html). For 3D examples, [LOBPCG](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html#rbbbc6164e7a5-4) is preferred since its convergence is faster.
3D examples are also considerably faster when using a GPU. GPU acceleration requires having [CuPy](https://docs.cupy.dev/en/stable/install.html) and [CUDA](https://developer.nvidia.com/cuda-downloads) installed in your computer. 

To use GPU acceleration in your 3D simulations, add the argument `method ='lobpcg-cupy'` in the Hamiltonian `solve` method. For example:

```python
eigenstates = H.solve( max_states = 50, method ='lobpcg-cupy')
```

The interface use [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) for input. In the file [constants.py](https://github.com/quantum-visualizations/qmsolve/blob/main/qmsolve/util/constants.py) there is a list of common conversion factors from other units, that can be imported and used to build your potential.

## Contributors

- [marl0ny](https://github.com/marl0ny)
- [Rafael de la Fuente](https://github.com/rafael-fuente)
- [Hudson Smith](https://github.com/dhudsmith)
