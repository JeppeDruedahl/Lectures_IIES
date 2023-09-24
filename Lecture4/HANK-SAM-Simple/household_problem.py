import numpy as np
import numba as nb

from consav.linear_interp import interp_1d_vec

@nb.njit
def solve_hh_backwards(par,z_trans,
    delta,lambda_u,w,RealR_ex_post,tau,
    vbeg_a_plus,vbeg_a,a,c,s,ss=False):
    """ solve backwards with vbeg_a_plus from previous iteration """
        
    # a. solution step
    for i_fix in range(par.Nfix):
        for i_z in range(par.Nz):
      
            # i. income
            if i_z == 0:
                y_pre = w
            else:
                y_pre = par.UI_ratio*w

            # ii. income after tax
            y = (1-tau)*y_pre

            # iii. EGM
            m = RealR_ex_post*par.a_grid + y

            # iv. consumption-saving
            if ss:

                c[i_fix,i_z,:] = 0.9*m
                a[i_fix,i_z,:] = m-c[i_fix,i_z,:]

            else:

                # o. EGM
                c_endo = (par.beta_grid[i_fix]*vbeg_a_plus[i_fix,i_z])**(-1/par.sigma)
                m_endo = c_endo + par.a_grid
            
                # oo. interpolation to fixed grid
                interp_1d_vec(m_endo,par.a_grid,m,a[i_fix,i_z])

                # ooo. enforce borrowing constraint
                a[i_fix,i_z,:] = np.fmax(a[i_fix,i_z,:],0.0)

                # oooo. implied consumption
                c[i_fix,i_z] = m-a[i_fix,i_z]

    # b. search (exogenous for now)
    for i_fix in range(par.Nfix):
        for i_z in range(par.Nz):
            if i_z == 0:
                s[i_fix,i_z,:] = delta
            else:
                s[i_fix,i_z,:] = 1.0

    fill_z_trans(par,z_trans,delta,lambda_u)

    # c. expectation step
    for i_fix in range(par.Nfix):
        for i_z_lag in range(par.Nz):
            for i_a_lag in range(par.Na):
                vbeg_a[i_fix,i_z_lag,i_a_lag] = 0.0    
                for i_z in range(par.Nz):
                    v_a = RealR_ex_post*c[i_fix,i_z,i_a_lag]**(-par.sigma)
                    vbeg_a[i_fix,i_z_lag,i_a_lag] += z_trans[i_fix,i_a_lag,i_z_lag,i_z]*v_a

################
# fill_z_trans #
################

@nb.njit(fastmath=True)
def fill_z_trans(par,z_trans,delta,lambda_u):
    """ transition matrix for z """
    
    for i_fix in nb.prange(par.Nfix):
        for i_a in nb.prange(par.Na):
            for i_z in nb.prange(par.Nz):
                for i_z_plus in nb.prange(par.Nz):

                    if i_z == 0:
                        
                        if i_z_plus == 0:
                            u_trans = 1.0-delta*(1-lambda_u)
                        else:
                            u_trans = delta*(1-lambda_u)
                    
                    else:

                        if i_z_plus == 0:
                            u_trans = lambda_u
                        else:
                            u_trans = 1.0-lambda_u

                    z_trans[i_fix,i_a,i_z,i_z_plus] = u_trans