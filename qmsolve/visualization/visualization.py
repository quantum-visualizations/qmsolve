from abc import abstractmethod 

# For type hinting purposes
from typing import Union, Dict, List


class Visualization:
    @abstractmethod
    def __init__(self,eigenstates):
        pass

    @abstractmethod
    def plot_eigenstate(self):
        pass

    @abstractmethod
    def slider_plot(self):
        pass
        
    @abstractmethod
    def animate_eigenstates(self):
        pass

    @abstractmethod
    def superpositions(self, states: Union[int, List[int], Dict[int, complex]],
                       **kw):
        """
        states - specify which eigenstates to superimpose on each other.
        If it's a single number then superimpose eigenstates 0 to states-1.
        kw - additional parameters to pass to the function.
        """
        pass


from .single_particle_1D import VisualizationSingleParticle1D
from .single_particle_2D import VisualizationSingleParticle2D
from .two_identical_particles_1D import VisualizationIdenticalParticles1D


def init_visualization(eigenstates):
    if eigenstates.type == "SingleParticle1D":
        return VisualizationSingleParticle1D(eigenstates) 

    elif eigenstates.type == "SingleParticle2D":
        return VisualizationSingleParticle2D(eigenstates) 

    elif eigenstates.type == "SingleParticle3D":
        from .single_particle_3D import VisualizationSingleParticle3D
        return VisualizationSingleParticle3D(eigenstates)

    elif eigenstates.type == "TwoIdenticalParticles1D":
        return VisualizationIdenticalParticles1D(eigenstates)
