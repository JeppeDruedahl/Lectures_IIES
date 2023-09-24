# find steady state

import time
import numpy as np
from scipy import optimize

from EconModel import jit

from consav.grids import equilogspace
from consav.markov import find_ergodic
from consav.misc import elapsed

import household_problem

def set_z_trans_ss(model):
    """ set z_trans """

    ss = model.ss

    with jit(model) as model_jit:
        household_problem.update_s(model_jit.par,ss.lambda_u_s,ss.s,ss.v_)
        household_problem.fill_z_trans(model_jit.par,ss.z_trans,ss.delta,ss.lambda_u_s,ss.s)

def get_Dz(model):
    """ get distribution for z """

    ss = model.ss

    return find_ergodic(ss.z_trans[0,0])

def set_Dbeg_ss(model):
    """ set initial distirubtion """

    par = model.par
    ss = model.ss

    Dz = get_Dz(model)
    for i_fix in range(par.Nfix):
        ss.Dbeg[i_fix,:,0] = par.beta_shares[i_fix]*Dz 
        ss.Dbeg[i_fix,:,1:] = 0.0      

def prepare_hh_ss(model):
    """ prepare the household block for finding the steady state """

    par = model.par
    ss = model.ss

    ############
    # 1. grids #
    ############

    # a. a
    par.a_grid[:] = equilogspace(0.0,par.a_max,par.Na)

    # b. z       
    par.z_grid[:] = 1.0
    
    #############################################
    # 2. transition matrix initial distribution #
    #############################################

    # a. transition matrix
    set_z_trans_ss(model)

    # b. ergodic
    set_Dbeg_ss(model)
    
    ################################################
    # 3. initial guess for intertemporal variables #
    ################################################

    model.set_hh_initial_guess()

def find_ss(model,do_print=False,calib_beta=True,recalib_beta=False):
    """ find the steady state"""

    t0 = time.time()

    assert not (calib_beta and recalib_beta), 'calib_beta and recalib_beta cannot be True at the same time'

    
    find_ss_SAM(model,do_print=do_print)
    find_ss_HANK(model,do_print=do_print,calib_beta=calib_beta,recalib_beta=recalib_beta)

    if do_print: print(f'steady state found in {elapsed(t0)}')

def find_ss_SAM(model,do_print=False):
    """ find the steady state - SAM """

    par = model.par
    ss = model.ss

    mean_UE = np.mean(model.data['UE'])/100
    mean_EU = np.mean(model.data['EU'])/100
    
    # a. shocks
    ss.shock_TFP = 1.0
    ss.retention_subsidy = 0.0
    ss.hiring_subsidy = 0.0

    # b. fixed
    ss.theta = par.theta_ss
    ss.w = par.w_ss
    ss.px = (par.epsilon_p-1)/par.epsilon_p
    ss.M = ss.px*ss.shock_TFP-ss.w

    ss.lambda_u = ss.lambda_u_s = mean_UE
    ss.delta = mean_EU

    # c. direct implications
    par.A = ss.lambda_u_s/ss.theta**(1-par.alpha)
    ss.lambda_v = par.A*ss.theta**(-par.alpha)

    # d. labor market dynamics
    ss.u = ss.delta/(ss.lambda_u_s+ss.delta)
    ss.vt = ss.u*ss.theta
    ss.v = (1-ss.lambda_v)*ss.vt
    ss.entry = ss.vt-(1-ss.delta)*ss.v

    # e. job and vacancy bellmans
    if par.exo_sep:

        par.p = np.nan
        par.Upsilon = np.nan
        ss.mu = 0.0

        ss.Vj = ss.M/(1-par.beta_firm*(1-ss.delta))

    else:

        par.p = par.p_fac*ss.delta

        Vj_Upsilon = (ss.delta/par.p)**(-1/par.psi)

        _nom = par.p*Vj_Upsilon**(-1)
        if np.abs(par.psi-1.0) < 1e-8:
            _nom *= np.log(Vj_Upsilon)
        else:
            _nom *= par.psi/(par.psi-1)*(1-Vj_Upsilon**(1-par.psi))

        _denom = (1-par.p*Vj_Upsilon**(-par.psi))

        mu_Vj = _nom/_denom

        ss.Vj = ss.M/(1+par.beta_firm*mu_Vj-par.beta_firm*(1-ss.delta))
        
        par.Upsilon = ss.Vj/Vj_Upsilon
        ss.mu = mu_Vj*ss.Vj

    if par.free_entry:
    
        par.kappa = ss.lambda_v*ss.Vj
        ss.Vv = 0.0
    
    else:

        _fac = 1-par.beta_firm*(1-ss.lambda_v)*(1-ss.delta)
        
        ss.Vv = par.kappa_0
        par.kappa = ss.lambda_v*ss.Vj - _fac*ss.Vv

    # f. misc
    ss.Y = ss.shock_TFP*(1-ss.u)
    ss.u_target = ss.u

def find_ss_HANK(model,do_print=False,calib_beta=False,recalib_beta=False):
    """ find the steady state - HANK """

    par = model.par
    ss = model.ss

    # a. shocks
    ss.shock_beta = 1.0

    # b. fixed
    ss.RealR_ex_post = ss.RealR = par.RealR_ss
    ss.phi_obar = par.phi_obar_ss
    ss.u_bar = par.u_bar_ss
    ss.tau = par.tau_ss
    ss.Pi = 1.0
    ss.R = ss.RealR*ss.Pi

    # c. dividends
    ss.adj_Vj = (1-ss.u)*ss.mu
    ss.adj_Vv = (ss.entry*ss.Vv)/(1+par.xi)
    ss.adj_Pi = par.phi/2*(ss.Pi-ss.Pi)**2*ss.Y
    ss.adj = ss.adj_Vj+ss.adj_Vv+ss.adj_Pi

    ss.div = ss.Y-(ss.w-ss.retention_subsidy)*(1-ss.u)-(1-par.adj_virtual_share)*ss.adj

    # d. household behavior
    ss.hh_div = 0.0 
    ss.hh_transfer = 0.0

    class StopCalibration(Exception): pass
    if not (calib_beta or recalib_beta):

        if par.exo_search:
            
            model.solve_hh_ss(do_print=do_print)
            model.simulate_hh_ss(do_print=do_print)

        else:

            def obj(x):

                par.zeta = x

                if do_print: print(f'{par.zeta = :.5f} -> ',end='')

                model.solve_hh_ss(do_print=False)
                model.simulate_hh_ss(do_print=False)

                if do_print: print(f'{ss.U_ALL_hh-ss.u = :5.1e}')

                return ss.U_ALL_hh-ss.u

            res = optimize.root_scalar(obj,bracket=[0.1,10.0],method='brentq')
            
            obj(res.root)            

        if do_print: print('')

    elif recalib_beta:

        assert par.Nfix == 3, 'not implemented'

        def obj(x):

            # i. unpack
            par.beta_max = par.beta_grid_min + x[0]*(par.beta_grid_max-par.beta_grid_min)

            # ii. implied beta grid and shares
            model.create_beta_grid()

            if do_print: print(f'x = [{x[0]:.5f},]: ',end='')

            # iii. moments
            model.solve_hh_ss()
            model.simulate_hh_ss()
            model.calc_moms_ss()

            # iv. error
            C_drop_ss = model.moms['C_drop_ss']
            error = (C_drop_ss-par.C_drop_ss_target)**2

            if do_print: print(f'{error:12.8f}, [{C_drop_ss = :.1f}]')
                
            if error < par.eps_calib_ss:
                raise StopCalibration 

            return error

        x0 = np.array([par.beta_max])
        bounds = ((0.0,1.0),)

        try:
            res = optimize.minimize(obj,x0,bounds=bounds,method='Nelder-Mead')
            objval = f'{obj(res.x):.1f}'
        except StopCalibration:
            if do_print: print('calibration stopped')
            objval = f'< {par.eps_calib_ss:.1f}'

    else:

        assert par.Nfix == 3, 'not implemented'

        # objective for calibration
        def obj(x):

            # i. unpack
            par.HtM_share = x[0]
            par.beta_max = par.beta_grid_min + x[1]*(par.beta_grid_max-par.beta_grid_min)

            # ii. implied beta grid and shares
            model.create_beta_grid()

            # iii. moments
            model.solve_hh_ss()
            model.simulate_hh_ss()
            model.calc_moms_ss()

            # iv. error
            error = 0.0

            C_drop_ss = model.moms['C_drop_ss']
            C_drop_ex_ss = model.moms['C_drop_ex']
            MPC_qtr = model.moms['MPC_qtr']

            error += (C_drop_ss-par.C_drop_ss_target)**2
            error += (C_drop_ex_ss-par.C_drop_ex_ss_target)**2

            if do_print:            
                print(f'x = [{x[0]:.5f},{x[1]:.5f}]: {error:12.8f} ',end='')
                print(f'[{MPC_qtr = :.1f}, {C_drop_ss = :.1f}, {C_drop_ex_ss = :.1f}]')

            if error < par.eps_calib_ss:
                raise StopCalibration 

            return error

        assert par.beta_max > par.beta_min
        
        x0_beta_max = (par.beta_max-par.beta_grid_min)/(par.beta_grid_max-par.beta_grid_min)
        x0 = np.array([par.HtM_share,x0_beta_max])
        bounds = ((0.0,1.0),(0.0,1.0))

        try:
            res = optimize.minimize(obj,x0,bounds=bounds,method='Nelder-Mead')
            objval = f'{obj(res.x):.1f}'
        except StopCalibration:
            if do_print: print('calibration stopped')
            objval = f'< {par.eps_calib_ss:.1f}'

    assert par.beta_max < par.beta_grid_max
    
    if do_print:

        if recalib_beta or calib_beta: print('')
        print(f'{par.HtM_share = : .3}')
        print(f'{par.PIH_share = : .3}')
        print(f'{1-par.HtM_share-par.PIH_share = : .3}')
        print(f'{par.beta_max**12 = :.3f}')
        print(f'{ss.U_UI_hh/ss.u*100 = :.2f} [{par.UI_of_u_share_target = :.2f}]')   
        if recalib_beta or calib_beta: print(f'calibration obj.: {objval}')
        print('')
        model.calc_moms_ss(vec=True,do_print=do_print)
        C_drop_ss = model.moms['C_drop_ss']
        C_drop_ex_ss = model.moms['C_drop_ex']
        MPC_qtr = model.moms['MPC_qtr']
        print(f'{MPC_qtr = :.1f}')
        print(f'{C_drop_ss = :.1f}')
        print(f'{C_drop_ex_ss = :.1f}')
        print('')

    # implied seachers and government bonds
    ss.S = ss.S_hh
    ss.qB = ss.A_hh

    # e. unemployment
    ss.U_UI_hh_guess = ss.U_UI_hh

    Dz = np.sum(ss.Dbeg,axis=2)
    u_ss = np.sum((par.i_u_hh > 0)*Dz)
    assert np.isclose(u_ss,ss.u)
    assert np.isclose(ss.U_ALL_hh,ss.u)
    assert ss.U_UI_hh <= ss.u

    if par.exo_search: assert np.isclose(ss.U_UI_hh/ss.u*100,par.UI_of_u_share_target)

    # f. pre-tax household income
    ss.UI = ss.phi_obar*ss.w*ss.U_UI_hh + par.phi_ubar*ss.w*(ss.u-ss.U_UI_hh)
    ss.Yt_hh = ss.w*(1-ss.u) + ss.UI
    ss.Y_hh = (1-ss.tau)*ss.Yt_hh

    # g. government budget
    ss.q = 1/(ss.RealR-par.delta_q)
    ss.B = ss.qB/ss.q

    expenses_no_G = ((1+par.delta_q*ss.q)*ss.B - ss.q*ss.B) + ss.UI
    taxes = ss.tau*ss.Yt_hh
    ss.G = taxes - expenses_no_G

    # h. clearing_Y
    ss.C_cap = ss.div
    C_tot = ss.C_hh + ss.C_cap + ss.G + (1-par.adj_virtual_share)*ss.adj
    ss.clearing_Y = ss.Y-C_tot

    if do_print:

        print(f'{ss.delta = :.3f}')
        print(f'{ss.lambda_u = :.3f}')
        print(f'{ss.qB/ss.Y_hh = :.3f}')
        print(f'{ss.qB/ss.Y = :.3f}')
        print(f'{ss.div/ss.Y = :.3f}')
        print(f'{ss.C_hh/ss.Y = :.3f}')
        print(f'{ss.C_cap/ss.Y = :.3f}')
        print(f'{ss.G/ss.Y = :.3f}')
        print(f'{ss.clearing_Y = :12.8f}')
