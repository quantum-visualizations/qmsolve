class Eigenstates:
    """[summary]
    """
    def __init__(self, energies, array, extent, N, type):
        """Info about the eigenstates

        Parameters
        ----------
        energies : array-like
            Energy eigenstates
        array : array-like
            The corresponding eigenvectors
        extent : array-like
            Extent of the simulation
        N : int
            Number of points to use for each dimension
        type : [type]
            [description] TODO
        """
        self.energies = energies
        self.array = array
        self.number = len(array)
        self.extent = extent
        self.N = N
        self.type = type
