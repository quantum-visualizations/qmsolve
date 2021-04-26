class Eigenstates:
    def __init__(self, energies, array, H, type):
        """Info about the eigenstates"""
        self.energies = energies
        self.array = array
        self.number = len(array)
        self.extent = H.extent
        self.N = H.N
        self.type = type
