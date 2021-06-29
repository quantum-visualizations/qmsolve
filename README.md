# Qmsolve: A module for solving and visualizing the Schrödinger equation.

This is an attempt at making a robust, easy to use solver, capable of solving and visualizing the Schrödinger equation for multiple particles, and representing the solutions both in 1D, 2D, and 3D.

The way this simulator works is by discretizing the Hamiltonian of an arbitrary potential and diagonalizing it for getting the energies and eigenstates of the system. It also features the possibility of interactively visualizing a superposition of the computed eigenstates and studying the time dependence of the resulting wavefunction.

This is a work in progress. Stay up to date about the next features!

## Installation

Just clone or download this repo.
The package requirements are:

1. numpy
2. matplotlib
3. scipy
4. mayavi (only for 3D simulations)

## Examples

Just run from the command line the corresponding Python scripts:

```
python 1D_harmonic_oscillator.py
```

![animation](/images/1D_harmonic_oscillator.gif)

```
python 3D_two_gaussian_wells.py
```

![animation](/images/3D_two_gaussian_wells.gif)

```
python 3D_four_gaussian_wells.py
```

![animation](/images/3D_four_gaussian_wells.gif)

```
python 1D_interactive_fermions_HO.py
```

![animation](/images/1D_interactive_fermions.gif)

```
python 1D_non_interactive_fermions_HO.py
```

![animation](/images/1D_non_interactive_fermions.gif)

In the examples from above you can check how in the non interactive case the energy levels are equally spaced and degenerated, however in the interactive case the degeneracy is broken.
As a starting point I suggest you to modify the confinement and the interaction potential to see what happens!

The interface use [Hartree atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) for input. In the file [constants.py](https://github.com/quantum-visualizations/qmsolve/blob/main/qmsolve/util/constants.py) there is a list of common conversion factors from other units, that can be imported and used to build your potential.

3D examples are considerably faster when using a GPU. GPU acceleration requires having [CuPy](https://docs.cupy.dev/en/stable/install.html) and [CUDA](https://developer.nvidia.com/cuda-downloads) installed in your computer. 

To use GPU acceleration in your 3D simulations, add the argument `method ='lobpcg-cupy'` in the Hamiltonian `solve` method. For example:

```python
eigenstates = H.solve( max_states = 50, method ='lobpcg-cupy')
```

## Contributors

- [marl0ny](https://github.com/marl0ny)
- [Rafael de la Fuente](https://github.com/rafael-fuente)
- [Hudson Smith](https://github.com/dhudsmith)
