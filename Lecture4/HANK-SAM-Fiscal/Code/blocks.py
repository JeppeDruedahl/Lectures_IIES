import numpy as np
import numba as nb

from GEModelTools import lag, lead, prev, next

@nb.njit
def delta_func(Vj,par,ss):    
    """ separations """

    if par.exo_sep:
        return ss.delta*np.ones(Vj.shape)
    else:
        return par.p*(np.fmax(Vj/par.Upsilon,1))**(-par.psi)

@nb.njit
def mu_func(Vj,par):
    """ continuation costs """

    if par.exo_sep:
        
        return np.zeros(Vj.shape)
    
    else:

        if np.abs(par.psi-1.0) < 1e-8:
            _fac = np.log(np.fmax(Vj/par.Upsilon,1.0))
        else:
            _fac = par.psi/(par.psi-1)*(1.0-(np.fmax(Vj/par.Upsilon,1.0))**(1-par.psi))

        _nom = par.p*par.Upsilon
        _denom = 1.0-par.p*(np.fmax(Vj/par.Upsilon,1.0))**(-par.psi)

        return _fac*_nom/_denom   

@nb.njit
def wage_and_profits(par,ini,ss,shock_TFP,px,w,retention_subsidy,M):

    if par.wage_setting == 'fixed':
        w[:] = ss.w

    M[:] = shock_TFP*px-(w-retention_subsidy)

@nb.njit
def Vj(par,ini,ss,Vj,M,delta,mu,errors_Vj):

    delta[:] = delta_func(Vj,par,ss)
    mu[:] = mu_func(Vj,par)

    Vj_plus = lead(Vj,ss.Vj)
    delta_plus = lead(delta,ss.delta)
    mu_plus = lead(mu,ss.mu)

    cont_Vj = (1-delta_plus)*par.beta_firm*Vj_plus - par.beta_firm*mu_plus
    errors_Vj[:] = Vj-(M+cont_Vj)

@nb.njit
def entry(par,ini,ss,Vv,entry):

    if not par.free_entry:   
        entry[:] = ss.entry*(np.fmax(Vv[:]/ss.Vv,0.0))**par.xi

@nb.njit
def labor_market(par,ini,ss,shock_TFP,u_target,vt,S,entry,delta,theta,lambda_v,lambda_u,lambda_u_s,u,v,errors_vt,errors_u,errors_u_target,Y):

    theta[:] = vt/S

    lambda_v[:] = par.A*theta**(-par.alpha)
    lambda_u_s[:] = par.A*theta**(1-par.alpha)

    v[:] = (1-lambda_v)*vt

    u_lag = lag(ini.u,u)
    v_lag = lag(ini.v,v)

    errors_vt[:] = vt - ((1-ss.delta)*v_lag + entry)
    errors_u[:] = u - (u_lag-S*lambda_u_s+delta*(1-u_lag))
    errors_u_target[:] = u-u_target
    lambda_u[:] = lambda_u_s*S/u_lag

    Y[:] = shock_TFP*(1-u)

@nb.njit
def Vv(par,ini,ss,Vj,Vv,hiring_subsidy,lambda_v,errors_Vv):
    
    if par.free_entry:

        LHS = np.zeros(Vv.shape)
        Vv_plus = np.zeros(Vv.shape)

    else:

        LHS = Vv
        Vv_plus = lead(Vv,ss.Vv)
    
    Vj_new_hire = Vj + hiring_subsidy
    RHS = -par.kappa + lambda_v*Vj_new_hire + (1-lambda_v)*(1-ss.delta)*par.beta_firm*Vv_plus

    errors_Vv[:] = LHS-RHS

@nb.njit
def wage_rule(par,ini,ss,u,shock_TFP,w,errors_WageRule):

    if par.wage_setting == 'fixed':

        pass

    else:

        for t in range(par.T):

            w_lag = prev(w,t,ini.w)
            
            curr_w = (u[t]/ss.u)**par.eta_u*( (1.0-u[t])/(1.0-ss.u))**par.eta_e*(shock_TFP[t])**par.eta_TFP
            RHS = ss.w*(w_lag/ss.w)**par.rho_w*curr_w**(1-par.rho_w)

            errors_WageRule[t] = w[t]-RHS

@nb.njit
def phillips(par,ini,ss,px,Pi,shock_TFP,u,errors_Pi):

    LHS = 1-par.epsilon_p + par.epsilon_p*px

    Pi_plus = lead(Pi,ss.Pi)        
    shock_TFP_plus = lead(shock_TFP,ss.shock_TFP)
    u_plus = lead(u,ss.u)

    RHS = par.phi*(Pi-ss.Pi)*Pi - par.beta_firm*par.phi*((Pi_plus-ss.Pi)*Pi_plus*(shock_TFP_plus*u_plus)/(shock_TFP*u))

    errors_Pi[:] = LHS-RHS

@nb.njit
def central_bank(par,ini,ss,R,Pi,RealR):

    for t in range(par.T):
        
        R_lag = prev(R,t,ini.R)

        R[t] = ss.R*(R_lag/ss.R)**(par.rho_R)*(Pi[t]/ss.Pi)**(par.delta_pi*(1-par.rho_R))
        
        if t < par.T-1:
            RealR[t] = R[t]/Pi[t+1]
        else:
            RealR[t] = R[t]/ss.Pi

@nb.njit
def dividends(par,ini,ss,u,vt,lambda_v,Y,entry,Vv,mu,Pi,adj_Vj,adj_Vv,adj_Pi,adj,w,retention_subsidy,hiring_subsidy,div):

    u_lag = lag(ini.u,u)
    adj_Vj[:] = (1-u_lag)*mu

    adj_Vv[:] = (entry*Vv)/(1+par.xi)
    adj_Pi[:] = par.phi/2*(Pi-ss.Pi)**2*Y

    adj[:] = adj_Vj+adj_Vv-adj_Pi

    div[:] = Y-(w-retention_subsidy)*(1-u) + hiring_subsidy*lambda_v*vt - (1-par.adj_virtual_share)*adj

@nb.njit
def arbitrage(par,ini,ss,RealR,q):

    for k in range(par.T):
        t = par.T-1-k
        q_plus = next(q,t,ss.q)
        q[t] = (1+par.delta_q*q_plus)/RealR[t]

@nb.njit
def government(par,ini,ss,UI,phi_obar,U_UI_hh_guess,w,u,Yt_hh,Y_hh,tau,B,G,hh_transfer,retention_subsidy,hiring_subsidy,lambda_v,vt,q,qB):

    UI[:] = phi_obar*w*U_UI_hh_guess + par.phi_ubar*w*(u-U_UI_hh_guess)
    Yt_hh[:] = w*(1-u) + UI
    Y_hh[:] = (1-tau)*Yt_hh
 
    for t in range(par.T):

        B_lag = prev(B,t,ini.B)
        tau[t] = ss.tau + par.omega*ss.q*(B_lag-ss.B)/ss.Yt_hh

        expenses = UI[t].copy()
        expenses += G[t]
        expenses += hh_transfer[t] 
        expenses += retention_subsidy[t]*(1-u[t])
        expenses += hiring_subsidy[t]*lambda_v[t]*vt[t]

        taxes = tau[t]*Yt_hh[t]

        B[t] = ((1+par.delta_q*q[t])*B_lag + expenses - taxes)/q[t]

    qB[:] = q*B

@nb.jit
def hh_vars(par,ini,ss,RealR,RealR_ex_post,q,div,hh_div,qB,C_cap):

    if par.div_return:
        
        hh_div[:] = 0.0
        hh_div_return = par.div_hh*(div-ss.div)
        C_cap[:] = div - hh_div_return

        RealR_ex_post[0] = (1+par.delta_q*q[0]+hh_div_return[0])*ini.B/ini.A_hh
        RealR_ex_post[1:] = RealR[:-1] + hh_div_return[1:]/qB[:-1]

    else:

        hh_div[:] = par.div_hh*(div-ss.div)
        C_cap[:] = div - hh_div[:]

        RealR_ex_post[0] = (1+par.delta_q*q[0])*ini.B/ini.A_hh
        RealR_ex_post[1:] = RealR[:-1]

@nb.njit
def market_clearing(par,ini,ss,qB,A_hh,errors_assets,U_UI_hh_guess,U_UI_hh,errors_U_UI,adj,C_RA,Y,C_hh,C_cap,G,clearing_Y,U_ALL_hh,u,errors_search,shock_beta,RealR):

    errors_U_UI[:] = U_UI_hh_guess-U_UI_hh
    errors_search[:] = U_ALL_hh - u
    
    if par.RA:

        C_RA = Y - (C_cap + G + (1-par.adj_virtual_share)*adj)
        C_plus = lead(C_RA,ss.C_hh)

        beta_RA_ss = 1/ss.RealR
        beta_RA = shock_beta*beta_RA_ss
        Euler_left = RealR*beta_RA*C_plus**(-par.sigma)
        Euler_right = C_RA**(-par.sigma)
        errors_assets[:] = Euler_left-Euler_right

    else:
        
        errors_assets[:] = qB-A_hh
    
    C_tot = C_hh + C_cap + G + (1-par.adj_virtual_share)*adj
    clearing_Y[:] = Y-C_tot

