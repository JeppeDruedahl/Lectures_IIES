# find steady state

import time
import numpy as np
from root_finding import brentq

from EconModel import jit

from consav.grids import equilogspace
from consav.markov import tauchen, find_ergodic
from consav.misc import elapsed

import household_problem

def prepare_hh_ss(model):
    """ prepare the household block for finding the steady state """

    par = model.par
    ss = model.ss

    ############
    # 1. grids #
    ############

    # a. beta
    par.beta_grid = np.array([par.beta_HtM,par.beta_mid,par.beta_PIH])
    
    par.beta_shares = np.array([par.HtM_share,0.0,par.PIH_share])
    par.beta_shares[1] = 1-np.sum(par.beta_shares)

    # b. assets
    par.a_grid[:] = equilogspace(0.0,par.a_max,par.Na)

    # c. income
    par.z_grid[:] = np.ones(par.Nz) # not used
    
    ###########################
    # 2. initial distribution #
    ###########################
    
    for i_fix in range(par.Nfix):
        ss.Dbeg[i_fix,:,0] = np.array([1-ss.u,ss.u])*par.beta_shares[i_fix]
        ss.Dbeg[i_fix,:,1:] = 0.0
    
    ################################################
    # 3. initial guess for intertemporal variables #
    ################################################

    model.set_hh_initial_guess()

def find_ss(model,do_print=False,fix_RealR=False):
    """ find the steady state"""

    t0 = time.time()
    
    find_ss_SAM(model,do_print=do_print)
    find_ss_HANK(model,do_print=do_print)

    if do_print: print(f'steady state found in {elapsed(t0)}')

def find_ss_SAM(model,do_print=False):
    """ find the steady state - SAM """

    par = model.par
    ss = model.ss
    
    # a. shocks
    ss.shock_TFP = 1.0

    # b. inflation and interest rates
    ss.Pi = 1.0
    ss.RealR_ex_post = ss.RealR = par.RealR_ss
    ss.R = ss.RealR*ss.Pi
    ss.q = 1/(ss.RealR-par.delta_q)

    # b. fixed
    ss.delta = par.delta_ss
    ss.lambda_u = par.lambda_u_ss
    ss.theta = par.theta_ss
    ss.px = (par.epsilon_p-1)/par.epsilon_p
    ss.w = par.w_share_ss*ss.px

    # c. direct implications
    par.A = ss.lambda_u/ss.theta**(1-par.alpha)
    ss.lambda_v = par.A*ss.theta**(-par.alpha)

    # d. labor market dynamics
    ss.u = ss.delta*(1-ss.lambda_u)/(ss.lambda_u+ss.delta*(1-ss.lambda_u))
    ss.ut = ss.u/(1-ss.lambda_u)
    ss.vt = ss.ut*ss.theta
    ss.v = (1-ss.lambda_v)*ss.vt
    ss.entry = ss.vt-(1-ss.delta)*ss.v
    ss.S = ss.vt/ss.theta

    # e. job and vacancy bellmans
    ss.Vj = (ss.px*ss.shock_TFP-ss.w)/(1-1/ss.RealR*(1-ss.delta))
    par.kappa = ss.lambda_v*ss.Vj

    if do_print:

        print(f'{par.A = :6.4f}')
        print(f'{par.kappa = :6.4f}')
        print(f'{ss.w = :6.4f}')
        print(f'{ss.delta = :6.4f}')
        print(f'{ss.lambda_u = :6.4f}')
        print(f'{ss.lambda_v = :6.4f}')
        print(f'{ss.theta = :6.4f}')
        print(f'{ss.u = :6.4f}')
        print(f'{ss.ut = :6.4f}')
        print(f'{ss.S = :6.4f}')


def find_ss_HANK(model,do_print=False):
    """ find the steady state - HANK """

    par = model.par
    ss = model.ss

    # a. shocks
    pass

    # c. policies
    ss.UI = par.UI_ratio*ss.w*ss.u
    ss.qB = par.qB_share_ss*ss.w
    ss.B = ss.qB/ss.q
    ss.Yt_hh = ss.w*(1-ss.u) + ss.UI
    ss.tau = ((1+par.delta_q*ss.q)*ss.B+ss.UI-ss.q*ss.B)/ss.Yt_hh

    # d. dividends  
    ss.div = ss.shock_TFP*(1-ss.u)*(1-ss.w/ss.shock_TFP)
    par.div_tax = par.div_tax_share_ss*ss.div
    ss.G = par.div_tax    
    ss.p_eq = (ss.div-par.div_tax)/(ss.RealR-1)

    # e. households
    model.solve_hh_ss(do_print=True)
    model.simulate_hh_ss(do_print=True)
    print(f'{ss.A_hh = }')

    A_hh_vec = np.array([np.sum(ss.a[i_fix]*ss.D[i_fix])/np.sum(ss.D[i_fix]) for i_fix in range(par.Nfix)])

    nom = ss.qB + par.mutual_fund_share*ss.p_eq - par.HtM_share*A_hh_vec[0] - (1-par.HtM_share)*A_hh_vec[1]
    denom = A_hh_vec[2]-A_hh_vec[1]
    par.PIH_share = nom/denom
    print(f'{par.PIH_share =}')

    model.solve_hh_ss(do_print=True)
    model.simulate_hh_ss(do_print=True)
    print(f'{ss.A_hh = }')

    ss.errors_assets = ss.qB + par.mutual_fund_share*ss.p_eq - ss.A_hh

    # f. clearing_Y
    C_cap = ss.div-par.div_tax
    C_tot = ss.C_hh + C_cap + ss.G
    ss.clearing_Y = ss.shock_TFP*(1-ss.u)-C_tot        