#!/usr/bin/python

from qutip import *
import numpy as np

def motion_problem(Omega, w, det, nbar, nMax):
    eta = 2 * np.pi / 578e-9 * np.sqrt(1.05e-34 / (2 * 171 * 1.67e-27 * (w * 1e6)))

    I_s = qeye(2)
    I_m = qeye(nMax)
    n_m = num(nMax)
    ikx_m = 1j * eta * (create(nMax) + destroy(nMax))

    H = (w * tensor(I_s, n_m) +
         Omega / 2 * tensor(sigmap(), (-ikx_m).expm()) + tensor(sigmam(), ikx_m.expm()) +
         det * tensor(fock_dm(2, 1), I_m))
    psi0 = tensor(fock_dm(2, 0), thermal_dm(nMax, nbar))

    return H, psi0

def solve_motion(p, tlist):
    H, psi0 = p
    return tlist, mesolve(H, psi0, tlist).states

def motion2_problem(Omega, w, det, nbar, nMax):
    eta = 2 * np.pi / 578e-9 * np.sqrt(1.05e-34 / (2 * 171 * 1.67e-27 * (w * 1e6)))

    I_s = qeye(2)
    I_m = qeye(nMax)
    n_m = num(nMax)
    ikx_m = 1j * eta * (create(nMax) + destroy(nMax))

    H = (w * tensor(n_m, I_s) +
         Omega / 2 * tensor((-ikx_m).expm(), sigmap()) + tensor(ikx_m.expm(), sigmam()) +
         det * tensor(I_m, fock_dm(2, 1)))
    psi0 = tensor(thermal_dm(nMax, nbar), fock_dm(2, 0))

    return H, psi0

def solve_motion2(p, tlist):
    H, psi0 = p
    return tlist, mesolve(H, psi0, tlist).states
