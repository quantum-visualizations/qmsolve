from abc import abstractmethod 


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


from .single_particle_1D import VisualizationSingleParticle1D
from .single_particle_2D import VisualizationSingleParticle2D


def init_visualization(eigenstates):
    if eigenstates.type == "SingleParticle1D":
        return VisualizationSingleParticle1D(eigenstates) 

    elif eigenstates.type == "SingleParticle2D":
        return VisualizationSingleParticle2D(eigenstates) 