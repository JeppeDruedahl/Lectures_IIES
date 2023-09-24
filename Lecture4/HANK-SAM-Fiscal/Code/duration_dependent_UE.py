import numpy as np
from scipy import optimize

import steady_state

def calibrate_s_eff_scale(model,do_print=False):
    """ calibrate s_eff so micro-u matches macro-u """

    par = model.par
    ss = model.ss

    steady_state.find_ss_SAM(model)

    # a. slope
    s_eff_ = par.s_eff.copy()
    def obj_u(x):

        par.s_eff[:] = s_eff_*x

        steady_state.set_z_trans_ss(model)
        Dz = steady_state.get_Dz(model)

        u = np.sum(Dz[model.par.i_u_hh[0] > 0])
        diff_u = 100*(u-ss.u)

        return diff_u
            
    res = optimize.root_scalar(obj_u,bracket=(0.25,2.0),method='bisect')
    obj_u(res.root)
    
    # check
    steady_state.set_z_trans_ss(model)
    Dz = steady_state.get_Dz(model)

    S = 0.0
    S += Dz[0]*par.s_eff[0]
    i_z = 1
    for i_UI in range(par.NUI):
        for i_u in range(1,par.Nu+1):
            S += Dz[i_z]*par.s_eff[i_u]
            i_z += 1

    assert np.isclose(S-ss.u,0.0)      
        
    if do_print: print(f' s_eff_scale = {res.root:.3f} [match ss.u]]')

def calibrate_UI_prob(model,do_print=False):
    """ calibrate UI prob to match share of u without UI """

    par = model.par

    steady_state.find_ss_SAM(model)

    def obj(x):

        par.UI_prob = x 

        steady_state.set_z_trans_ss(model)
        Dz = steady_state.get_Dz(model)
        
        u = np.sum(Dz[model.par.i_u_hh[0] > 0])

        UI_of_u_share = 0.0
        for i_z in range(par.Nz):
                
            i_u = par.i_u_hh[0,i_z]
            i_UI = par.i_UI_hh[0,i_z]

            if i_UI == 1 and i_u >= 1 and i_u <= par.u_bar_ss:
                UI_of_u_share += Dz[i_z]

        UI_of_u_share /= u
        
        diff = par.UI_of_u_share_target-UI_of_u_share*100

        return diff

    res = optimize.root_scalar(obj,bracket=(0.1,0.9),method='bisect')
    obj(res.root)

    if do_print: print(f' {par.UI_prob = :.3f} [math par.UI_of_u_share_target]')    

def get_scale_and_UI_prob(model,do_print=False):
    """ calibrate s_eff_scale and UI_prob"""

    par = model.par

    it = 0
    while True:
        
        UI_prob = model.par.UI_prob
        if do_print: print(f'iteration {it}, {par.UI_prob = :.3f}')   

        calibrate_s_eff_scale(model,do_print=do_print)
        calibrate_UI_prob(model,do_print=do_print)

        diff = np.abs(UI_prob - model.par.UI_prob)
        if do_print: print(f'diff(UI_prob) = {diff:.2e}\n')   

        if it > 100 or diff < 1e-8:
            if do_print: print('converged') 
            break   

        it += 1

def calibrate(model,do_print=False):
    """ calibrate duration dependence of UE """

    # assumption 1: same u dynamics for all beta types
    # assumption 2: same duration dynamics for all UI outcomes

    par = model.par

    # a. initial guess
    par.UI_prob = 0.5
    par.s_eff[0] = 0.0
    par.s_eff[1:par.UEs_taget.size+1] = np.exp(par.UEs_taget)
    par.s_eff[par.UEs_taget.size+1:] = par.s_eff[par.UEs_taget.size]

    # b. calibrate
    get_scale_and_UI_prob(model,do_print=do_print)