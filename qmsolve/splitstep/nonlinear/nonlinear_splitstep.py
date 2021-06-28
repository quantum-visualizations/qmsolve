from .. import SplitStepMethod
import numpy as np
from typing import Union, Callable, Tuple
from ...util import constants


class NonlinearSplitStepMethod(SplitStepMethod):

    """
    Split-Operator method for the non-linear Schrodinger equation.
    
    References:

      Xavier Antoine, Weizhu Bao, Christophe Besse
      Computational methods for the dynamics of 
      the nonlinear Schrodinger/Gross-Pitaevskii equations.
      Comput. Phys. Commun., Vol. 184, pp. 2621-2633, 2013.
      https://arxiv.org/pdf/1305.1093

    """

    def __init__(self, hamiltonian, timestep):
        SplitStepMethod.__init__(self, hamiltonian, timestep)
        self._nonlinear = lambda psi: psi

    def __call__(self, psi: np.ndarray) -> np.ndarray:
        """
        Step the wavefunction in time.
        """
        psi = self._nonlinear(psi)
        psi_p = np.fft.fftn(psi*self._exp_potential)
        psi_p = psi_p*self._exp_kinetic
        psi = np.fft.ifftn(psi_p)*self._exp_potential
        psi = self._nonlinear(psi)
        if self._norm:
            psi = psi/np.sqrt(np.sum(psi*np.conj(psi)))
        return psi
    
    def set_nonlinear_term(self, nonlinear_func: Callable) -> None:
        """
        Set the nonlinear term.
        """
        self._nonlinear = nonlinear_func


class CoupledTwoSystemNonlinearSplitStepMethod(SplitStepMethod):
    
    def __init__(self, hamiltonian, potential, timestep, **kw):
        params = {'m1': constants.m_e, 'm2': constants.m_e,
                  'lambda1': 0.0, 'lambda2': 0.0,
                  'V1': potential, 'V2': potential,
                  'nonlinear1': lambda psi: psi,
                  'nonlinear2': lambda psi: psi,
                  'hbar': constants.hbar}
        for key in kw.keys():
            if key in params:
                params[key] = kw[key]
            else:
                raise KeyError
        self._m1 = params['m1']
        self._m2 = params['m2']
        self._lambda1 = params['lambda2']
        self._lambda2 = params['lambda1']
        self._V1 = params['V1']
        self._V2 = params['V2']
        self._nonlinear1 = params['nonlinear1']
        self._nonlinear2 = params['nonlinear2']
        self._hbar = params['hbar']
        self._exp_potential1 = None
        self._exp_potential2 = None
        self._exp_kinetic1 = None
        self._exp_kinetic2 = None
        self._exp_c = None
        SplitStepMethod.__init__(self, hamiltonian, timestep)
    
    def set_potential(self, V: np.ndarray, V2: np.ndarray = None) -> None:
        self._V1 = V
        self._V2 = V if not V2 else V2
    
    def set_nonlinear_term(self, nonlinear: Callable, 
                           nonlinear2: Callable = None) -> None:
        """
        Set the nonlinear term.
        """
        self._nonlinear1 = nonlinear
        self._nonlinear2 = nonlinear if not nonlinear2 else nonlinear2
    
    def set_timestep(self, timestep: Union[float, np.complex128]) -> None:
        """
        Set the timestep. It can be real or complex.

        The coupling term is derived using the following lines of code:

        >>> from sympy import Matrix, Symbol, exp
        >>> p2 = Symbol('p2', real=True)
        >>> lambda1 = Symbol('lambda1', real=True)
        >>> lambda2 = Symbol('lambda2', real=True)
        >>> hbar, dt = Symbol('hbar', positive=True), Symbol('dt')
        >>> ec = exp(-1.0j*dt*Matrix([[0, lambda1], 
                     [lambda2, 0]])/hbar).simplify()
        >>> for i, term in enumerate(ec):
        >>>     j, k = i // 2, i % 2
        >>>     print('e%d%d = ' % (j, k), term)


        """
        self._dt = timestep
        dt, hbar = self._dt, self._hbar
        lambda1, lambda2 = self._lambda1, self._lambda2
        m1, m2 = self._m1, self._m2
        dt_inv_hbar = self._dt/hbar
        I = 1.0j
        cosh, sinh = np.cosh, np.sinh
        self._exp_potential1 = np.exp(-0.25j*dt_inv_hbar*self._V1)
        self._exp_potential2 = np.exp(-0.25j*dt_inv_hbar*self._V2)
        p = np.meshgrid(*[2.0*np.pi*hbar*np.fft.fftfreq(d)*d/
                          self._dim[i] for i, d in enumerate(self.V.shape)])
        p2 = sum([p_i**2 for p_i in p])
        self._exp_kinetic1 = np.exp(-0.5j*(self._dt*p2/(2.0*m1*hbar)))
        self._exp_kinetic2 = np.exp(-0.5j*(self._dt*p2/(2.0*m2*hbar)))
        if lambda1 == 0.0 and lambda2 == 0.0:
            e00, e01, e10, e11 = 1.0, 0.0, 0.0, 1.0
            self._exp_c = [[e00, e01], [e10, e11]]
            return
        e00 = 1.0*cosh(dt*I*(lambda1*lambda2)**0.5/hbar)
        e01 = -1.0*lambda1*(lambda1*lambda2)**(-0.5)* \
            sinh(dt*I*(lambda1*lambda2)**0.5/hbar)
        e10 = -1.0*(lambda1*lambda2)**0.5* \
            sinh(dt*I*(lambda1*lambda2)**0.5/hbar)/lambda1
        e11 = 1.0*cosh(dt*I*(lambda1*lambda2)**0.5/hbar)
        self._exp_c = [[e00, e01], [e10, e11]]


    def __call__(self, psi1: np.ndarray, 
                 psi2: np.ndarray) -> Tuple[np.ndarray]:
        """
        Step the wavefunction in time.
        """
        psi1 = self._nonlinear1(psi1)*self._exp_potential1
        psi2 = self._nonlinear2(psi2)*self._exp_potential2
        psi1_p = np.fft.fftn(psi1)
        psi2_p = np.fft.fftn(psi2)
        psi1 = np.fft.ifftn(psi1_p*self._exp_kinetic1)
        psi2 = np.fft.ifftn(psi2_p*self._exp_kinetic2)
        psi1 = self._nonlinear1(psi1)*self._exp_potential1
        psi2 = self._nonlinear2(psi2)*self._exp_potential2
        psi1 = self._exp_c[0][0]*psi1 + self._exp_c[0][1]*psi2
        psi2 = self._exp_c[1][0]*psi1 + self._exp_c[1][1]*psi2
        if self._norm:
            psi1 = psi1/np.sqrt(np.sum(psi1*np.conj(psi1)))
            psi2 = psi2/np.sqrt(np.sum(psi2*np.conj(psi2)))
        return psi1, psi2
    
