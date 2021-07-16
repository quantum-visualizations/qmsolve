# QMsolve: A module for solving and visualizing the Schrödinger equation

![animation](/images/3D_two_gaussian_wells.gif)

QMsolve seeks to provide a solid and easy to use solver, capable of solving the Schrödinger equation for one and two particles, 
and creating descriptive and stunning visualizations of its solutions both in 1D, 2D, and 3D.


## Installation

```
pip install qmsolve
```

## How it works

The way this simulator works is by discretizing the Hamiltonian with an arbitrary potential, 
specified as a function of the particle observables. This is achieved with the `Hamiltonian` constructor.

Then, the `Hamiltonian.solve` method efficiently diagonalizes the Hamiltonian and outputs the energies and the eigenstates of the system.
Finally, the eigenstates can be plotted with the use of the `visualization` class.

The `visualization.superpositions` method features the possibility of interactively visualizing a superposition of the computed eigenstates and studying the time dependence of the resulting wavefunction. 

For a quick start, take a look to the examples found in the [examples subdirectory](https://github.com/quantum-visualizations/qmsolve/tree/main/examples).

## Examples

To perform the simulations, just run from the corresponding Python scripts on the command prompt.

```
python 1D_harmonic_oscillator.py
```

![animation](/images/1D_harmonic_oscillator.gif)

```
python 1D_interactive_fermions_HO.py
```

![animation](/images/1D_interactive_fermions.gif)

```
python 1D_non_interactive_fermions_HO.py
```

![animation](/images/1D_non_interactive_fermions.gif)

```
python 3D_four_gaussian_wells.py
```

![animation](/images/3D_four_gaussian_wells.gif)


In the two interactive particle examples from above you can check how in the non interactive case the energy levels are equally spaced and degenerated, however in the interactive case the degeneracy is broken.
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
