#!/usr/bin/env python

from qutip import *
import numpy as np
# from scipy.linalg import expm
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from qutip.superoperator import spre
from qutip.expect import expect_rho_vec
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import expm, expm_multiply

def _add_H(M, N, H):
    if isinstance(H, (QobjEvo,)):
        raise Exception("Time dependent Hamiltonian not supported.")
    H = Qobj(H)
    for i in range(N):
        for j in range(N):
            mi = i * N + j
            for k in range(N):
                M[mi, k * N + j] += H[i, k] / 1j
                M[mi, i * N + k] -= H[k, j] / 1j

def _add_C(M, N, C):
    if isinstance(C, (QobjEvo,)):
        raise Exception("Time dependent collapse term not supported.")
    C = Qobj(C)
    Ch = C.conj().trans()
    L = Ch * C
    for i in range(N):
        for j in range(N):
            mi = i * N + j
            for k in range(N):
                for l in range(N):
                    M[mi, k * N + l] += C[i, k] * Ch[l, j]
                M[mi, k * N + j] -= L[i, k] / 2
                M[mi, i * N + k] -= L[k, j] / 2

def _construct_matrix(H, c_ops, use_sparse):
    hshape = H.shape
    assert len(hshape) == 2
    assert hshape[0] == hshape[1]
    N = hshape[0]

    coherent = len(c_ops) == 0

    if coherent:
        if use_sparse:
            return csc_matrix(Qobj(H), dtype=complex), N, True
        else:
            return np.array(Qobj(H), dtype=complex), N, True

    if use_sparse:
        M = lil_matrix((N**2, N**2), dtype=complex)
    else:
        M = np.zeros((N**2, N**2), dtype=complex)
    _add_H(M, N, H)
    for C in c_ops:
        _add_C(M, N, C)
    if use_sparse:
        M = csc_matrix(M)
    return M, N, False

def const_mesolve(H, rho0, tlist, c_ops=None, e_ops=None, progress_bar=None,
                  use_sparse=True):
    if c_ops is None:
        c_ops = []
    elif isinstance(c_ops, (Qobj,)):
        c_ops = [c_ops]
    elif isinstance(c_ops, (QobjEvo,)):
        raise Exception("Time dependent collapse term not supported.")

    if e_ops is None:
        e_ops = []
    if isinstance(e_ops, Qobj):
        e_ops = [e_ops]

    if isinstance(e_ops, dict):
        e_ops_dict = e_ops
        e_ops = [e for e in e_ops.values()]
    else:
        e_ops_dict = None

    if progress_bar is None:
        progress_bar = BaseProgressBar()
    if progress_bar is True:
        progress_bar = TextProgressBar()

    if isket(rho0):
        rho0 = ket2dm(rho0)

    n_tsteps = len(tlist)
    output = solver.Result()
    output.solver = "const_mesolve"
    output.times = tlist
    output.expect = []
    output.states = []
    output.num_collapse = len(c_ops)

    e_ops_data = []
    if callable(e_ops):
        n_expt_op = 0
        expt_callback = True
        output.num_expect = 1
    elif isinstance(e_ops, list):
        n_expt_op = len(e_ops)
        expt_callback = False
        output.num_expect = n_expt_op
        for op in e_ops:
            if not isinstance(op, Qobj) and callable(op):
                output.expect.append(np.zeros(n_tsteps, dtype=complex))
                e_ops_data.append(None)
                continue
            if op.dims != rho0.dims:
                raise TypeError(f"e_ops dims ({op.dims}) are not "
                                f"compatible with the state's "
                                f"({rho0.dims})")
            e_ops_data.append(spre(op).data)
            if op.isherm and rho0.isherm:
                output.expect.append(np.zeros(n_tsteps))
            else:
                output.expect.append(np.zeros(n_tsteps, dtype=complex))
    else:
        raise TypeError("Expectation parameter must be a list or a function")

    M, N, coherent = _construct_matrix(H, c_ops, use_sparse)

    progress_bar.start(n_tsteps)

    if not coherent:
        rhof0 = np.ndarray.flatten(np.array(rho0))
        for t_idx, t in enumerate(tlist):
            progress_bar.update(t_idx)
            # TODO: try use the form of expm_multiply with sampling
            rhof = expm_multiply(t * M, rhof0)
            rho = Qobj(rhof.reshape(N, N))
            output.states.append(rho)

            if expt_callback:
                # use callback method
                output.expect.append(e_ops(t, rho))

            for m in range(n_expt_op):
                op = e_ops[m]
                if not isinstance(op, Qobj) and callable(op):
                    output.expect[m][t_idx] = op(t, rho_t)
                    continue
                output.expect[m][t_idx] = expect_rho_vec(e_ops_data[m], rhof,
                                                         op.isherm and rho0.isherm)
    else:
        _rho0 = np.array(rho0, dtype=complex)
        for t_idx, t in enumerate(tlist):
            progress_bar.update(t_idx)
            U = expm((t / 1j) * M)
            _rho = U @ _rho0 @ U.conj().T
            rho = Qobj(_rho)
            output.states.append(rho)
            rhof = np.ndarray.flatten(_rho)

            if expt_callback:
                # use callback method
                output.expect.append(e_ops(t, rho))

            for m in range(n_expt_op):
                op = e_ops[m]
                if not isinstance(op, Qobj) and callable(op):
                    output.expect[m][t_idx] = op(t, rho_t)
                    continue
                output.expect[m][t_idx] = expect_rho_vec(e_ops_data[m], rhof,
                                                         op.isherm and rho0.isherm)

    if e_ops_dict:
        output.expect = {e: output.expect[n] for n, e in enumerate(e_ops_dict.keys())}

    return output
