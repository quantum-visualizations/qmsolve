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


class TimeVisualization:
    @abstractmethod
    def __init__(self,simulation):
        pass

    @abstractmethod
    def plot(self):
        pass

    @abstractmethod
    def animate(self):
        pass



from .single_particle_1D import VisualizationSingleParticle1D
from .single_particle_2D import VisualizationSingleParticle2D
from .two_identical_particles_1D import VisualizationIdenticalParticles1D
from ..eigenstates import Eigenstates


def init_eigenstate_visualization(eigenstates):

    if eigenstates.type == "SingleParticle1D":
        return VisualizationSingleParticle1D(eigenstates) 

    elif eigenstates.type == "SingleParticle2D":
        return VisualizationSingleParticle2D(eigenstates) 

    elif eigenstates.type == "SingleParticle3D":
        from .single_particle_3D import VisualizationSingleParticle3D
        return VisualizationSingleParticle3D(eigenstates)

    elif eigenstates.type == "TwoIdenticalParticles1D":
        return VisualizationIdenticalParticles1D(eigenstates)

from .single_particle_1D import TimeVisualizationSingleParticle1D
from .single_particle_2D import TimeVisualizationSingleParticle2D
from .two_identical_particles_1D import TimeVisualizationTwoIdenticalParticles1D
from ..particle_system import SingleParticle, TwoParticles
from ..time_dependent_solver import TimeSimulation

def init_timesimulation_visualization(simulation):

    if (isinstance(simulation.H.particle_system ,SingleParticle) and (simulation.H.spatial_ndim == 1)):
        return TimeVisualizationSingleParticle1D(simulation) 

    elif (isinstance(simulation.H.particle_system , SingleParticle) and (simulation.H.spatial_ndim == 2)):
        return TimeVisualizationSingleParticle2D(simulation) 

    elif (isinstance(simulation.H.particle_system , SingleParticle) and (simulation.H.spatial_ndim == 3)):
        raise NotImplementedError()
        #from .single_particle_3D import TimeVisualizationSingleParticle3D
        #return TimeVisualizationSingleParticle3D(simulation)

    elif isinstance(simulation.H.particle_system , TwoParticles):
        return TimeVisualizationTwoIdenticalParticles1D(simulation)


def init_visualization(argument):

    if isinstance(argument, Eigenstates):
        eigenstates = argument
        return init_eigenstate_visualization(eigenstates)

    elif isinstance(argument, TimeSimulation):
        simulation = argument
        return init_timesimulation_visualization(simulation)