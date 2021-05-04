from abc import ABC, abstractmethod


class ParticleSystem(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_observables(self, H):
        pass

    @abstractmethod
    def get_kinetic_matrix(self, H):
        pass
        
    @abstractmethod
    def get_eigenstates(self, H, max_states, eigenvalues, eigenvectors):
        pass
