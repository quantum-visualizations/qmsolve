"""
Single particle quantum mechanics simulation
using the split-operator method.

References:
https://www.algorithm-archive.org/contents/
split-operator_method/split-operator_method.html

https://en.wikipedia.org/wiki/Split-step_method

"""
import numpy as np
from ..util.constants import hbar, Å, femtoseconds
from ..hamiltonian import Hamiltonian
import time
import matplotlib.pyplot as plt
from ..util.colour_functions import complex_to_rgba

from .split_step import SplitStep, SplitStepCupy


class TimeSimulation:
    """
    Class for configuring time dependent simulations.
    """

    def __init__(self, hamiltonian, method = "split-step"):

        self.H = hamiltonian

        implemented_solvers = ('split-step', 'split-step-cupy')

        if method == "split-step":
            self.method = SplitStep(self)
        elif method == "split-step-cupy":
            self.method = SplitStepCupy(self)
        else:
            raise NotImplementedError(
                f"{method} solver has not been implemented. Use one of {implemented_solvers}")



    def run(self, initial_wavefunction, total_time, dt, store_steps = 1):
        """
        """
        self.method.run(initial_wavefunction, total_time, dt, store_steps)