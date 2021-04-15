import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .two_particles import TwoIdenticalParticles
from ..util.constants import *


class TwoBosons(TwoIdenticalParticles):

    def __init__(self, m = m_e, spin = None):
        TwoIdenticalParticles.__init__(self, m=m_e, spin=None, p=1)
