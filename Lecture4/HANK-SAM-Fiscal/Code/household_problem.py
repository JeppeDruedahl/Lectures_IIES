import numpy as np
import numba as nb

from consav.linear_interp import interp_1d_vec

@nb.jit
def isclose(a,b):
    return np.abs(a-b) < 1e-12

@nb.njit
def solve_hh_backwards(par,z_trans,
    delta,lambda_u_s,shock_beta,w,RealR_ex_post,tau,u_bar,phi_obar,hh_div,hh_transfer,
    vbeg_plus,vbeg_a_plus,vbeg,vbeg_a,a,c,v_,s,u_ALL,u_UI,ss=False):
    """ solve backwards with vbeg_a_plus from previous iteration """

    # a. solution step
    vend = np.zeros((par.Nfix,par.Nz,par.Na)) # end-of-period value
    for i_fix in range(par.Nfix):
        for i_z in range(par.Nz):
            
            z = par.z_grid[i_z] # productivity
            i_u = par.i_u_hh[i_fix,i_z] # unemployment indicator
            i_UI = par.i_UI_hh[i_fix,i_z] # UI indicator

            # i. income
            if i_u == 0:
                u_UI_ = 0.0
                yt = w*z
                u_ALL[i_fix,i_z,:] = 0.0
            elif i_UI == 0:
                u_UI_ = 0.0
                yt = par.phi_ubar*w*z
                u_ALL[i_fix,i_z,:] = 1.0
            else:
                u_UI_ = np.fmax(np.fmin(u_bar-(i_u-1),1.0),0.0)
                yt = (u_UI_*phi_obar + (1-u_UI_)*par.phi_ubar)*w*z
                u_ALL[i_fix,i_z,:] = 1.0

            u_UI[i_fix,i_z,:] = u_UI_

            # ii. income after tax
            if par.Nfix == 3 and par.PIH_share > 0.0:                
                if i_fix == 2:
                    y = (1-tau)*yt + par.div_PIH*hh_div/par.PIH_share + hh_transfer
                else:
                    y = (1-tau)*yt + (1-par.div_PIH)*hh_div/(1-par.PIH_share) + hh_transfer
            else:
                y = (1-tau)*yt + hh_div + hh_transfer

            # iii. EGM
            m = RealR_ex_post*par.a_grid + y

            # iv. consumption-saving
            if par.beta_grid[i_fix] < par.beta_HtM_cutoff or isclose(par.beta_shares[i_fix],0.0): # HtM
        
                a[i_fix,i_z,:] = 0.0
                c[i_fix,i_z,:] = m
                if par.do_vbeg: v_[i_fix,i_z,:] = c[i_fix,i_z]**(1-par.sigma)/(1-par.sigma)

            elif ss:

                c[i_fix,i_z,:] = 0.9*m
                a[i_fix,i_z,:] = m-c[i_fix,i_z,:]
                if par.do_vbeg: v_[i_fix,i_z,:] = c[i_fix,i_z]**(1-par.sigma)/(1-par.sigma)

            else:

                # o. EGM
                if par.do_vbeg: vend_endo = shock_beta*par.beta_grid[i_fix]*vbeg_plus[i_fix,i_z]
                c_endo = (shock_beta*par.beta_grid[i_fix]*vbeg_a_plus[i_fix,i_z])**(-1/par.sigma)
                m_endo = c_endo + par.a_grid

                # oo. interpolation to fixed grid
                interp_1d_vec(m_endo,par.a_grid,m,a[i_fix,i_z])
                if par.do_vbeg: interp_1d_vec(m_endo,vend_endo,m,vend[i_fix,i_z])

                # ooo. enforce borrowing constraint
                a[i_fix,i_z,:] = np.fmax(a[i_fix,i_z,:],0.0)

                # oooo. implied consumption and value
                c[i_fix,i_z] = m-a[i_fix,i_z]
                if par.do_vbeg: 
                    v_[i_fix,i_z] = c[i_fix,i_z]**(1-par.sigma)/(1-par.sigma) + vend[i_fix,i_z]
                    if i_z == 0: v_[i_fix,i_z] -= par.vartheta

    # b. update transition matrix
    update_s(par,lambda_u_s,s,v_)
    fill_z_trans(par,z_trans,delta,lambda_u_s,s)

    # c. expectation step
    v_a = RealR_ex_post*c**(-par.sigma)
    for i_fix in range(par.Nfix):
        for i_z_lag in range(par.Nz):
            
            i_u_lag = par.i_u_hh[i_fix,i_z_lag]
            s_eff = par.s_eff[np.int_(i_u_lag)]            
            vbeg[i_fix,i_z_lag,:] = 0.0
            if par.do_vbeg and i_u_lag > 0:
                vbeg[i_fix,i_z_lag,:] -= par.zeta*(s[i_fix,i_z,:]/s_eff)**(1+par.nu)/(1+par.nu)

            vbeg_a[i_fix,i_z_lag,:] = 0.0
            
            for i_z in range(par.Nz):

                if par.do_vbeg: vbeg[i_fix,i_z_lag,:] += z_trans[i_fix,:,i_z_lag,i_z]*v_[i_fix,i_z,:]
                vbeg_a[i_fix,i_z_lag,:] += z_trans[i_fix,:,i_z_lag,i_z]*v_a[i_fix,i_z,:]

################
# fill_z_trans #
################

@nb.njit
def update_s(par,lambda_u_s,s,v):

    for i_fix in range(par.Nfix):

        s[i_fix,0,:] = 0.0

        i_z = 0
        for i_UI in range(par.NUI):
            for i_z in range(par.Nz):

                # a. search efficiency
                i_u = par.i_u_hh[i_fix,i_z]
                s_eff = par.s_eff[np.int_(i_u)]

                # b. exogenous search
                if par.exo_search or isclose(par.beta_shares[i_fix],0.0):

                    s[i_fix,i_z] = s_eff

                # c. endogenous search
                else:
                    
                    for i_a in range(par.Na):
                        
                        i_u_u = np.fmin(np.int_(i_u)+1,par.Nu-1)
                        i_z_u = 1 + i_UI*par.Nu + i_u_u
                        
                        v_u = v[i_fix,np.int_(i_z_u),i_a]
                        v_e = v[i_fix,0,i_a]
                        v_diff = v_e-v_u

                        s_ = (lambda_u_s*v_diff/par.zeta)**(1/par.nu)
                        s_ = np.fmax(s_,0.01)
                        s_ = np.fmin(s_,10.0)

                        s[i_fix,i_z,i_a] = s_eff*s_

                i_z += 1

@nb.njit
def fill_z_trans(par,z_trans,delta,lambda_u_s,s):
    """ transition matrix for z """
    
    # b. transition matrix
    for i_fix in range(par.Nfix):
        for i_a_lag in range(par.Na):
            for i_z_lag in range(par.Nz):
                for i_z in range(par.Nz):
                            
                    i_u_lag = par.i_u_hh[i_fix,i_z_lag] 
                    i_UI_lag = par.i_UI_hh[i_fix,i_z_lag]

                    i_u = par.i_u_hh[i_fix,i_z] 
                    i_UI = par.i_UI_hh[i_fix,i_z]                            

                    if i_u_lag == 0: # working last period
                        
                        if i_u == 0:
                            z_trans_ = 1.0-delta
                        elif i_u == 1 and i_UI == 0:
                            z_trans_ = delta*(1-par.UI_prob)
                        elif i_u == 1 and i_UI == 1:
                            z_trans_ = delta*par.UI_prob
                        else:
                            z_trans_ = 0.0
                    
                    else: # unemployed last

                        lambda_u_s_ = lambda_u_s*s[i_fix,i_z_lag,i_a_lag]

                        if i_u == 0:
                            z_trans_ = lambda_u_s_
                        elif not i_UI_lag == i_UI:
                            z_trans_ = 0.0
                        elif (i_u == i_u_lag+1) or (i_u == i_u_lag == par.Nu):
                            z_trans_ = 1.0-lambda_u_s_
                        else:
                            z_trans_ = 0.0

                    z_trans_ = np.fmin(z_trans_,1.0)
                    z_trans_ = np.fmax(z_trans_,0.0)

                    z_trans[i_fix,i_a_lag,i_z_lag,i_z] = z_trans_