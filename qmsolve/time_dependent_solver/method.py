from abc import ABC, abstractmethod


class Method(ABC):
    @abstractmethod
    def __init__(self,simulation):
        pass

    @abstractmethod
    def run(self, H):
        pass
