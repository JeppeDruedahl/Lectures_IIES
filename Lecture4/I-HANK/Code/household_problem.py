import numpy as np
import numba as nb

from consav.linear_interp import interp_1d_vec
from GEModelTools import isclose

import grids

@nb.njit(parallel=True, cache=False)
def solve_hh_backwards_1A(par, ra, UniformT, tau, wnT, wnNT, eps_beta, d, AI_lag, z_trans, 
                          vbeg_a_plus, vbeg_a, a, c, c_T, c_NT, inc_T,inc_NT, uc_T, uc_NT):
    """ solve backwards with vbeg_a from previous iteration - one asset"""

    # a. EGM loop
    for i_beta in nb.prange(par.Nbeta):
        for i_z in nb.prange(par.Nz):

            if par.s_set[i_z] == 0:
                w = wnT
                T = UniformT * par.sT_vec[0]                
            elif par.s_set[i_z] == 1:
                w = wnNT
                T = UniformT * par.sT_vec[1]

            m = (1 + ra) * par.a_grid + (1-tau) * w * par.z_grid_ss[i_z] + T + d

            # EGM
            c_endo = (eps_beta * (par.beta_grid[i_beta]) * vbeg_a_plus[i_beta, i_z]) ** (-1 / par.CRRA)
            m_endo = c_endo + par.a_grid

            # interpolation
            interp_1d_vec(m_endo, par.a_grid, m, a[i_beta, i_z, :])

            # enforce borrowing constraint
            a[i_beta, i_z, :] = np.fmax(a[i_beta, i_z, :], par.a_min)
            c[i_beta, i_z, :] = m - a[i_beta, i_z, :]

        # expectation step
        va = (1 + ra) * c[i_beta] ** (-par.CRRA)
        vbeg_a[i_beta] = z_trans[i_beta] @ va

    # b. subtract from C premium
    c[:] = c - par.xi*AI_lag

    # c. extra output
    c_T[:, :par.Ne, :] = c[:, :par.Ne, :] 
    c_NT[:, par.Ne:, :] = c[:, par.Ne:, :] 

    for i_z in nb.prange(par.Nz):
        if par.s_set[i_z] == 0:
            inc_T[:, i_z, :] = a[:, i_z, :] * ra + (1-tau)*wnT + UniformT * par.sT_vec[0]
            uc_T[:, i_z, :] = c[:,i_z, :]**(-par.CRRA) * par.z_grid_ss[i_z]
        elif par.s_set[i_z] == 1:
            inc_NT[:, i_z, :] = a[:, i_z, :] * ra + (1-tau)*wnNT + UniformT * par.sT_vec[1]
            uc_NT[:, i_z, :] = c[:,i_z, :]**(-par.CRRA) * par.z_grid_ss[i_z]

@nb.njit(parallel=True, cache=False)
def solve_hh_backwards_2A(par, ra, rai, UniformT, tau, wnT, wnNT, eps_beta, z_trans, vbeg_a_plus, 
                          vbeg_a, a, ai, d_payout, c, c_T, c_NT, inc_T,inc_NT, uc_T, uc_NT):
    """ solve backwards with vbeg_a from previous iteration - two assets"""

    # unpack
    AI_target = par.AI_target
    rai_ss = par.rai_ss_target

    # a. EGM loop
    for i_beta in nb.prange(par.Nbeta):
        for i_z in nb.prange(par.Nz):

            if par.s_set[i_z] == 0:
                labor_inc = wnT
                T = UniformT * par.sT_vec[0]
            elif par.s_set[i_z] == 1:
                labor_inc = wnNT
                T = UniformT * par.sT_vec[1]
            
            # loop over illiquid asset 
            for i_ai_lag in range(par.Nai): # end-of-previous-period

                # EGM
                c_endo = (eps_beta * (par.beta_grid[i_beta]) * vbeg_a_plus[i_beta, i_z,:,i_ai_lag]) ** (-1 / par.CRRA)
                m_endo = c_endo + par.a_grid

                d_payout[i_beta, i_z, :, i_ai_lag] = rai_ss/(1+rai_ss)*par.ai_grid[i_ai_lag]  \
                                            + par.chi*(par.ai_grid[i_ai_lag] - (1+rai_ss)*AI_target)
                
                income = (1-tau)*labor_inc * par.z_grid_ss[i_z] + T
                m = (1 + ra) * par.a_grid + income  + d_payout[i_beta, i_z, :, i_ai_lag]

                # liquid assets
                interp_1d_vec(m_endo, par.a_grid, m, a[i_beta, i_z, :, i_ai_lag])

                # enforce borrowing constraint
                a[i_beta, i_z, :, i_ai_lag] = np.fmax(a[i_beta, i_z, :, i_ai_lag], par.a_min)
                c[i_beta, i_z, :, i_ai_lag] = m - a[i_beta, i_z, :, i_ai_lag]

                # illiquid assets
                ai[i_beta, i_z, :, i_ai_lag] = (1+rai) * par.ai_grid[i_ai_lag] - d_payout[i_beta, i_z, :, i_ai_lag]

        
        # expectation step
        for i_ai_lag in range(par.Nai):  
            va = (1 + ra) * c[i_beta,:,:,i_ai_lag] ** (-par.CRRA)
            vbeg_a[i_beta,:,:,i_ai_lag] = z_trans[i_beta] @ va

    # b. subtract illiquid premium from C
    for i_beta in nb.prange(par.Nbeta):
        for i_z in nb.prange(par.Nz):
            for i_ai_lag in range(par.Nai): # end-of-previous-period
                c[i_beta, i_z, :, i_ai_lag] = c[i_beta, i_z, :, i_ai_lag] - par.xi*par.ai_grid[i_ai_lag]

    # c. extra output
    c_T[:, :par.Ne, :,:] = c[:, :par.Ne, :,:]
    c_NT[:, par.Ne:, :,:] = c[:, par.Ne:, :,:]

    for i_z in nb.prange(par.Nz):
        if par.s_set[i_z] == 0:
            inc_T[:, i_z, :,:] = a[:, i_z, :,:] * ra + (1-tau)*wnT + UniformT * par.sT_vec[0]
            uc_T[:, i_z, :] = c[:,i_z, :]**(-par.CRRA) * par.z_grid_ss[i_z]
        elif par.s_set[i_z] == 1:
            inc_NT[:, i_z, :,:] = a[:, i_z, :,:] * ra + (1-tau)*wnNT + UniformT * par.sT_vec[1]
            uc_NT[:, i_z, :] = c[:,i_z, :]**(-par.CRRA) * par.z_grid_ss[i_z]

@nb.njit
def util(par, c, N, psi):
    """ utility function """

    if isclose(par.sigma, 1.0):
        u = np.log(c) - psi * (N) ** (1 + par.inv_frisch) / (1 + par.inv_frisch)
    else:
        u = (c) ** (1 - par.CRRA) / (1 - par.CRRA) - psi * (N) ** (1 + par.inv_frisch) / (1 + par.inv_frisch)

    return u

def compute_RA_jacs(ss, par):
    """ compute jacobians for RA model """

    U = np.triu(np.ones((par.T, par.T)), k=0)

    beta = par.beta_mean
    for j in range(par.T):
        par.M_Y[j, :] = (1 - beta) * beta ** np.arange(par.T)

    par.M_R[:] = -1 / par.CRRA * (np.eye(par.T) - par.M_Y) @ U * ss.CR * beta
    par.M_beta[:] = -1 / par.CRRA * (np.eye(par.T) - par.M_Y) @ U * (1 + ss.r) * ss.CR

    # # Impose steady state in final period
    # # Note: Should be done differently if we want to allow for permament shocks to C 
    par.M_R[-1,:] = 0. 
    par.M_beta[-1,:] = 0. 
    par.M_Y[-1,:] = 0. 

    if par.HH_type == 'TA-IM':

        par.M_R[:] *= (1 - par.MPC_est[0])
        par.M_Y[:] *= (1 - par.MPC_est[0])
        par.M_Y[:] += par.MPC_est[0] * np.eye(par.T)
    

def prepare_hh_ss(model):
    """ prepare the household block for finding the steady state """

    par = model.par
    ss = model.ss

    ############
    # 1. grids #
    ############

    grids.create_grids_full(model)

    #############################################
    # 2. transition matrix initial distribution #
    #############################################

    if par.HH_type == 'HA':
        for i_beta in range(par.Nbeta):
            ss.Dbeg[i_beta, :, 0] = par.z_ergodic_ss / par.Nbeta 
            ss.Dbeg[i_beta, :, 1:] = 0.0
    else:
        ss.Dbeg[:] = par.z_ergodic_ss[np.newaxis,:,np.newaxis,np.newaxis]
        #ss.Dbeg[i_beta, :, 1:,1:] = 0.0     
        ss.Dbeg[:] = ss.Dbeg / np.sum(ss.Dbeg.flatten())       

    ################################################
    # 3. initial guess for intertemporal variables #
    ################################################

    if par.HH_type == 'HA':
        vbeg_a = np.zeros((par.Nfix, par.Nz, par.Na))
        m = (1 + ss.ra) * par.a_grid[np.newaxis, :] + par.z_grid_ss[:, np.newaxis] + ss.UniformT
    else:
        vbeg_a = np.zeros((par.Nfix, par.Nz, par.Na, par.Nai))
        m = (1 + ss.ra) * par.a_grid[np.newaxis, :, np.newaxis] + \
            par.z_grid_ss[:, np.newaxis, np.newaxis] + ss.UniformT + 0.05*par.ai_grid[np.newaxis, np.newaxis,:]

    a = 0.90 * m  # pure guess
    c = m - a
    vbeg_a[:] = (1 + ss.ra) * c[np.newaxis, ...] ** (-par.CRRA)

    for i_fix in range(par.Nfix):
        ss.vbeg_a[i_fix] = vbeg_a[i_fix]
