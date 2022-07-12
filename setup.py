from setuptools import setup, find_packages

long_description = """
# QMsolve: A module for solving and visualizing the Schrödinger equation
QMsolve seeks to provide easy solid and easy to use solver, capable of solving the Schrödinger equation for one and two particles, 
and creating descriptive and stunning visualizations of its solutions both in 1D, 2D, and 3D.

## Installation

```
pip install qmsolve
```

## How it works

The way this simulator works is by discretizing and Hamiltonian of an arbitrary potential, 
specified as a function of the particle observables. This is achieved with the `Hamiltonian`constructor.
Then, the `Hamiltonian.solve` method efficiently diagonalizes the Hamiltonian and outputs the energies and the eigenstates of the system.
Finally, the eigenstates can be plotted with the `visualization` class.

The `visualization.superpositions` method features the possibility of interactively visualizing a superposition of the computed eigenstates and studying 
the time dependence of the resulting wavefunction. 

For efficiently solving the time dependent Schrödinger equation, `TimeSimulation` class must be used. It takes as an argument the Hamiltonian you have previously defined, and the method you desire to use.

For a quick start, take a look to the examples found in the examples subdirectory.
"""


setup(
    name='qmsolve',
    version='2.1.0',
    description='A module for solving and visualizing the Schrödinger equation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/quantum-visualizations/qmsolve',
    download_url='https://github.com/quantum-visualizations/qmsolve/archive/main.zip',
    keywords = ['schrödinger-equation', 'quantum-physics', 'quantum-mechanics', 'quantum-programming'],
    author='Rafael de la Fuente, marl0ny, Hudson Smith',
    author_email='rafael.fuente.herrezuelo@gmail.com, lammarlonaurel123456789@gmail.com, dhudsmith@gmail.com',
    license='BSD-3-Clause',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'matplotlib', 'h5py', 'progressbar'],
    extras_require = {'with_mayavi:python_version=="3.7"': ['vtk == 8.1.2', 'mayavi >= 4.7.1', 'PyQt5 >= 5.14.1'],
                      'with_mayavi:python_version>="3.8"': ['vtk == 9.0.1', 'mayavi >= 4.7.1', 'PyQt5 >= 5.14.1'] },
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',

    ],
    include_package_data = True,
    python_requires ='>=3.7',
)