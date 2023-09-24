
import numpy as np
import numba as nb

from GEModelTools import lag, lead, isclose

###########
# helpers #
###########

@nb.njit
def annual_inflation(pi_ann, P, T, ssP):
    for t in range(T):
        if t >=4:
            pi_ann[t] = P[t]/P[t-4] - 1
        else:
            pi_ann[t] = P[t]/ssP - 1

@nb.njit
def convert_q2a(x_q):
    x_a = (1+x_q)**4 - 1  
    return x_a 

@nb.njit
def I_adj_cost(I_p, I, phi_I):
    adj_cost = phi_I/2 * (I_p/I-1)**2
    adj_cost_deriv = phi_I * (I_p/I-1)
    return adj_cost, adj_cost_deriv

@nb.njit
def CTW_adj_cost(I_p, I, Spp):
    x = I_p/I
    sqrt_Spp = np.sqrt(Spp)
    adj_cost = 1/2 * (np.exp(sqrt_Spp*(x-1)) + np.exp(-sqrt_Spp*(x-1)) - 2) 
    adj_cost_deriv = sqrt_Spp/2 * (np.exp(sqrt_Spp*(x-1)) - np.exp(-sqrt_Spp*(x-1))) 
    return adj_cost, adj_cost_deriv

@nb.njit
def chain_price_index(P, ssP, Y, P1, ssP1, Y1, P2, ssP2, Y2, T, sign=1.0):
    for t in range(T):
        if t==0:
            Y[t]  = (ssP1 * Y1[t] + sign*(ssP2 * Y2[t]))/ssP
        else:
            Y[t]  = (P1[t-1] * Y1[t] + sign*(P2[t-1] * Y2[t]))/P[t-1]
        P[t] = (P1[t] * Y1[t] + sign*(P2[t] * Y2[t]))/Y[t]

@nb.njit
def CES_demand(c, PH, PF, P, eta, alpha, c_bar=0.0):
    cF = alpha     * (PF/P)**(-eta) * (c - c_bar*PF/P) + c_bar  
    cH = (1-alpha) * (PH/P)**(-eta) * (c - c_bar*PF/P)
    return cF, cH

@nb.njit
def CES_demand_3_inputs(X, PH, PNT, PF,PX, eta, Theta_T, Theta_NT, Theta_M):
    XT2T  = Theta_T     * (PH/PX)**(-eta) * X 
    XNT2T = Theta_NT * (PNT/PX)**(-eta) * X
    XM2T  = Theta_M * (PF/PX)**(-eta) * X
    return XT2T, XNT2T, XM2T

@nb.njit
def get_intermediates(XT, XNT, PXT, PXNT, PH, PNT, PF, Theta_T, Theta_NT, etaX):  
    XT2T, XNT2T, XM2T = CES_demand_3_inputs(XT, PH, PNT, PF, PXT, etaX, Theta_T[0], Theta_T[1], Theta_T[2]) # intermediate demand in tradeable sector
    XNT2NT, XT2NT = CES_demand(XNT, PH, PNT, PXNT, etaX, Theta_NT[0]) # intermediate demand in non-tradeable sector
    return XT2NT, XNT2T, XNT2NT, XT2T, XM2T

@nb.njit
def CES_demand_T(C, PT, PNT, P, etaT, alphaT, c_bar=0.0):
    CT = alphaT     * (PT/P)**(-etaT) * (C - c_bar*PT/P) + c_bar  
    CNT = (1-alphaT) * (PNT/P)**(-etaT) * (C - c_bar*PNT/P) 
    return CT, CNT

@nb.njit
def Armington(PH_s, PF_s, C_s, gamma, alpha):
    CH_s = alpha * (PH_s/PF_s)**(-gamma) * C_s
    return CH_s

@nb.njit
def Armington_dyn(PH_s, PF_s, C_s, gamma, alpha, phi_X, T, PH_s_ss, PF_s_ss):
    
    ToT = PH_s/PF_s
    ToT_tilde = np.zeros_like(ToT)
    ToT_tilde[0] = ToT[0]**(1-phi_X)*(PH_s_ss/PF_s_ss)**phi_X
    for t in range(1,T):
        ToT_tilde[t] = ToT[t]**(1-phi_X)*ToT_tilde[t-1]**phi_X
    
    CH_s = alpha*(ToT_tilde)**(-gamma) * C_s

    return CH_s

@nb.njit
def Price_index(PH, PF, eta, alpha):
    if isclose(eta,1.0):
        P = PF**alpha * PH**(1-alpha)
    else:
        P = (alpha*PF**(1-eta) + (1-alpha)*PH**(1-eta))**(1/(1-eta))
    return P

@nb.njit
def Price_index_T(PH, PNT, PF, etaX, X_share_T):
    if isclose(etaX,1.0):
        PXT = PH**X_share_T[0] * PNT**X_share_T[1] * PF**X_share_T[2]
    else:
        PXT = (X_share_T[0]*PH**(1-etaX) + X_share_T[1]*PNT**(1-etaX) + X_share_T[2]*PF**(1-etaX))**(1/(1-etaX))
    return PXT

@nb.njit
def Price_index_NT(PH, PNT, etaX, X_share_NT):
    if isclose(etaX,1.0):
        PXNT = PNT**X_share_NT[0] * PH**X_share_NT[1]  
    else:
        PXNT = (X_share_NT[0]*PNT**(1-etaX) + X_share_NT[1]*PH**(1-etaX))**(1/(1-etaX))
    return PXNT

@nb.njit
def sol_Price_index_rel(pF, eta, alpha):
    if isclose(eta,1.0):
        pH = (1.0/(pF**alpha))**(1/(1-alpha))
    else:
        pH = ((1.0-alpha*pF**(1-eta))/(1-alpha))**(1/(1-eta))
    return pH

@nb.njit
def sol_Price_index2(pH, pF, eta, alpha):
    if isclose(eta,1.0):
        pF = (1.0/(pH**(1-alpha)))**(1/alpha)
    else:
        pF = ((1.0-(1-alpha)*pH**(1-eta))/alpha)**(1/(1-eta))
    return pF

@nb.njit
def Inf(P, ssval):
    P_lag = lag(ssval,P) 
    pi = P/P_lag - 1
    pi_p = lead(pi,0.)
    return pi, pi_p

@nb.njit
def P_from_inf(P, inf, T, ssP):
    for t in range(T):
        if t == 0:
            P[t] = ssP*(1+inf[t]) 
        else:
            P[t] = P[t-1]*(1+inf[t]) 

@nb.njit
def Get_HH_A(A, C, I, ra, T, ss):
    for t in range(T):  
        if t==0:
            A[t] = ss.A * (1+ra[t]) + I[t]   - C[t] 
        else:
            A[t] = A[t-1] * (1+ra[t]) + I[t]   - C[t]    

@nb.njit    
def sol_backwards_lin(var, ssvar, a, b, T):
    for tt in range(T):
        t = T-1-tt
        if t == T-1:
            var[t] = ssvar
        else:
            var[t] = a[t] + b[t] * var[t+1]   
            

###############
#    blocks   #
###############

# block for variables not used in transition paths 
@nb.njit
def misc_block(par, ini, ss, PTP, PO_T, PO_NT, OT, ONT, goods_mkt, AI, Atot, AI_lag, Y_s): 
    pass

@nb.njit
def risk_premium_block(RP_T, RP_NT, mT, mNT, par, ss):
    RP_T[:] = 1 - par.mu_r * np.tanh(par.kappa_r*(mT / ss.mT - 1))
    RP_T_p = lead(RP_T,ss.RP_T)
    RP_T_pp = lead(RP_T_p,ss.RP_T)

    RP_NT[:] = 1 - par.mu_r * np.tanh(par.kappa_r*(mNT / ss.mNT - 1))
    RP_NT_p = lead(RP_NT,ss.RP_NT)
    RP_NT_pp = lead(RP_NT_p,ss.RP_NT)

    return RP_T_p, RP_T_pp, RP_NT_p, RP_NT_pp


@nb.njit
def investment_block_I_adj(par, ss, KT, KNT, IT, INT, I, rKT, rKNT, QT, 
                            QNT, r, IT_FOC, INT_FOC, KTp, KNTp,
                            adjI_T, adjI_NT, ZT, ZNT, PH, PNT, mcT, 
                            mcNT, rHT, rHNT, HT, HNT, P, VAT, RP_T, RP_NT,
                            mT, mNT, ITp, INTp, subI):

    #adj_cost_f = CTW_adj_cost
    adj_cost_f = I_adj_cost

    if par.use_capital:
        # Get risk premium 
        RP_T_p, RP_T_pp, RP_NT_p, RP_NT_pp = risk_premium_block(RP_T, RP_NT, mT, mNT, par, ss)

        # Get investment
        IT[:] = lag(ss.ITp, ITp)
        INT[:] = lag(ss.INTp, INTp)
        adjI_T[:] = adj_cost_f(IT, lag(ss.IT, IT), par.phi_I)[0] * IT
        adjI_NT[:] = adj_cost_f(INT, lag(ss.INT, INT), par.phi_I)[0] * INT

        # Use law of motion to obtain capital from investment
        for t in range(par.T):
            if t == 0:
                KT[t] = ss.KT
                KNT[t] = ss.KNT
            else:
                KT[t] = (1-par.delta_K) * KT[t-1] + IT[t-1]
                KNT[t] = (1-par.delta_K) * KNT[t-1] + INT[t-1]

        # Get return on capital from capital demand
        rKT[:] = rHT * (KT/(HT * (par.TFP * par.TFP_s[0]) ** (par.sigma_NK - 1))  / (par.prodalphaK[0] ** (par.sigma_NK)))**(-1/par.sigma_NK)    
        rKNT[:] = rHNT * (KNT/(HNT * (par.TFP * par.TFP_s[1]) ** (par.sigma_NK - 1))  / (par.prodalphaK[1] ** (par.sigma_NK)))**(-1/par.sigma_NK)   

        # Obtain Tobins Q through backwards iteration
        QT[-1] = 1.
        QNT[-1] = 1.

        for tt in range(par.T-1):
            t = par.T-1-tt-1
            QT[t] = (rKT[t+1] * RP_T_p[t]     + (1-par.delta_K) * QT[t+1]) / (1+r[t])
            QNT[t] = (rKNT[t+1] * RP_NT_p[t]   + (1-par.delta_K) * QNT[t+1]) / (1+r[t])

        # generate various leads
        r_p = lead(r, ss.r)
        QT_p = lead(QT, ss.QT)
        QNT_p = lead(QNT, ss.QNT)
        subI_p = lead(subI,ss.subI)

        ITp_lag = IT
        INTp_lag = INT
        ITpp = lead(ITp, ss.IT)
        INTpp = lead(INTp, ss.INT)

        # Investment focs
        adj_T,d_adj_T = adj_cost_f(ITp, ITp_lag, par.phi_I)
        adj_NT, d_adj_NT = adj_cost_f(INTp, INTp_lag, par.phi_I)
        d_adj_T_p = adj_cost_f(ITpp, ITp, par.phi_I)[1]
        d_adj_NT_p = adj_cost_f(INTpp, INTp, par.phi_I)[1]

        IT_FOC[:] = (1-subI_p +  adj_T + d_adj_T*(ITp/ITp_lag)) - (QT_p/RP_T_p + RP_T_pp/RP_T_p * 1/(1+r_p) * (d_adj_T_p * (INTpp/INTp)**2))
        INT_FOC[:] = (1-subI_p + adj_NT + d_adj_NT*(INTp/INTp_lag)) - (QNT_p/RP_NT_p + RP_NT_pp/RP_NT_p * 1/(1+r_p) * (d_adj_NT_p * (INTpp/INTp)**2))

        # Aggregate investment
        I[:] = lead(IT + INT, ss.I) # Should be leaded since I is predetermined in model, but not in data/LP
    else:
        rKT[:] = 0.
        rKNT[:] = 0.
        KT[:] = 0.
        KNT[:] = 0.
        IT[:] = 0.
        INT[:] = 0.
        ITp[:] = 0.
        INTp[:] = 0.
        I[:] = 0.
        IT_FOC[:] = 0.
        INT_FOC[:] = 0.
        adjI_T[:] = 0.
        adjI_NT[:] = 0.

@nb.njit
def get_labor_adj_cost(r, ss_r, N, ss_N, phi_N):

    N_l = lag(ss_N, N)
    N_p = lead(N,ss_N)

    adj_term = phi_N*(N/N_l-1) + 1/(1+ss_r)*(phi_N/2 * (N_p/N-1)**2 - phi_N * (N_p/N-1)*N_p/N)

    return adj_term


@nb.njit
def production_block(par, ini, ss, ZT, ZNT, Z, mcT, mcNT, XT, XNT, XNT2T, XT2NT, XNT2NT, 
                    XT2T, XM2T, PXT, PXNT, subP,subI, VAT, piH, piNT, PH, PNT, PF, rHT, rHNT, HT, HNT, 
                    WT, WNT, rWT, rWNT, NT, NNT, N, PGDP, GDP, YT, YNT, Y, wnT,
                    wnNT, w, P, wT, wNT, KT, KNT, rKT, rKNT, r,
                    IT, INT, I, QT, QNT, IT_FOC, INT_FOC, KTp, KNTp, RP_T, RP_NT, 
                    mT, mNT,  ITp, INTp, adjI_T, adjI_NT, labor_comp, Lfoc_T, Lfoc_NT):

    PC_beta = 1 / (1 + ss.r)

    # Solve NKPCs     
    ZT_p = lead(ZT, ss.ZT)
    ZNT_p = lead(ZNT, ss.ZNT)
    VAT_p = lead(VAT, VAT[-1])
    piH_p = lead(piH, ss.piH)
    piNT_p = lead(piNT, ss.piNT)
    piH_l = lag(ss.piH, piH)
    piNT_l = lag(ss.piNT, piNT)

    # invert NKPC to get real marginal costs
    disc_T = PC_beta * ZT_p / ZT * (1 - VAT_p) / (1 - VAT)
    inf_term_T = ((1+piH)/(1+piH_l)**par.pi_index - 1) * ((1+piH)/(1+piH_l)**par.pi_index)  
    inf_term_p_T = ((1+piH_p)/(1+piH)**par.pi_index - 1) * ((1+piH_p)/(1+piH)**par.pi_index)  
    mcT[:]  = 1/par.mu[0]  + (inf_term_T - disc_T * inf_term_p_T)/(par.NKslope[0]*PH/P)  

    disc_NT = PC_beta * ZNT_p / ZNT * (1 - VAT_p) / (1 - VAT)
    inf_term_NT = ((1+piNT)/(1+piNT_l)**par.pi_index - 1) * ((1+piNT)/(1+piNT_l)**par.pi_index)  
    inf_term_p_NT = ((1+piNT_p)/(1+piNT)**par.pi_index - 1) * ((1+piNT_p)/(1+piNT)**par.pi_index)  
    mcNT[:]  = 1/par.mu[1]  + (inf_term_NT - disc_NT * inf_term_p_NT)/(par.NKslope[1]*PNT/P)  
  

    # CES prices 
    PXT[:] = Price_index_T(PNT, PH, PF, par.etaX, par.X_share_T)
    PXNT[:] = Price_index_NT(PH, PNT, par.etaX, par.X_share_NT)

    if isclose(par.sigma_XH, 1.0) and ss.XT > 0:  # Cobb-Douglas
        rHT[:]  = (mcT * PH * (1 - VAT) / (
                    (PXT * (1 - VAT) / par.prodalphaX[0]) ** par.prodalphaX[0])) ** (1 / (1 - par.prodalphaX[0])) * (
                            1 - par.prodalphaX[0])  / P
        rHNT[:] = (mcNT * PNT * (1 - VAT) / (
                    (PXNT * (1 - VAT) / par.prodalphaX[1]) ** par.prodalphaX[1])) ** (1 / (1 - par.prodalphaX[1])) * (
                            1 - par.prodalphaX[1])  / P
    elif isclose(ss.XT, 0.0):  # Linear production function with labor
        rHT[:] = (1 - VAT) * mcT * PH   / P # check TFP here
        rHNT[:] = (1 - VAT) * mcNT * PNT / P # check TFP here
    else:  # CES
        rHT[:] = ((((1 - VAT) * mcT * PH * par.TFP * par.TFP_s[0]) ** (1 - par.sigma_XH) - par.prodalphaX[
            0] ** par.sigma_XH * PXT ** (1 - par.sigma_XH)) / par.prodbeta[0] ** par.sigma_XH) ** (
                            1 / (1 - par.sigma_XH))  / P
        rHNT[:] = ((((1 - VAT) * mcNT * PNT * par.TFP * par.TFP_s[1]) ** (1 - par.sigma_XH) - par.prodalphaX[
            1] ** par.sigma_XH * PXNT ** (1 - par.sigma_XH)) / par.prodbeta[1] ** par.sigma_XH) ** (
                             1 / (1 - par.sigma_XH)) / P

    # Intermediate goods
    if par.alphaX[0] > 0.0 and par.alphaX[1] > 0.0:
        XT[:] = par.prodalphaX[0] ** (par.sigma_XH) * (PXT / (PH * mcT)) ** (-par.sigma_XH) * ZT 
        XNT[:] = par.prodalphaX[1] ** (par.sigma_XH) * (PXNT / (PNT * mcNT)) ** (-par.sigma_XH) * ZNT 
        XT2NT[:], XNT2T[:], XNT2NT[:], XT2T[:], XM2T[:] = get_intermediates(XT, XNT, PXT, PXNT, PH, PNT, PF,
                                                                                    par.X_share_T, par.X_share_NT,                                                                      par.etaX)
    else:
        XT[:] = 0.0
        XNT[:] = 0.0
        XNT2T[:] = 0.0
        XT2NT[:] = 0.0
        XT2NT[:], XNT2T[:], XNT2NT[:], XT2T[:] = 0.0, 0.0, 0.0, 0.0
        PXT[:], PXNT[:] = 0.0, 0.0
        XM2T[:] = 0.

    # Next stage in production function 
    # Invert production function to get H given Z and X 
    exp = (par.sigma_XH - 1) / par.sigma_XH
    if isclose(par.sigma_XH, 1.0) and ss.XT > 0:
        HT[:] = (ZT / (XT ** par.alphaX[0])) ** (1 / (1 - par.alphaX[0]))
        HNT[:] = (ZNT / (XNT ** par.alphaX[1])) ** (1 / (1 - par.alphaX[1]))
    elif ss.XT > 0:
        HT[:] = ((ZT  ** exp - par.prodalphaX[0] * XT ** exp) / par.prodbeta[0]) ** (1 / exp)
        HNT[:] = ((ZNT ** exp - par.prodalphaX[1] * XNT ** exp) / par.prodbeta[1]) ** (1 / exp)
    else:
        HT[:] = ZT
        HNT[:] = ZNT

    # Given rH and H get solve for investment and capital 
    investment_block_I_adj(par, ss, KT, KNT, IT, INT, I, rKT, rKNT, QT, 
                            QNT, r, IT_FOC, INT_FOC, KTp, KNTp,
                            adjI_T, adjI_NT, ZT, ZNT, PH, PNT, mcT, 
                            mcNT, rHT, rHNT, HT, HNT, P, VAT, RP_T, RP_NT, 
                            mT, mNT, ITp, INTp, subI)

    # Invert production function for H to get N
    exp_NK = (par.sigma_NK - 1) / par.sigma_NK
    if par.use_capital:
        if isclose(par.sigma_NK, 1.0):
            NT[:] = (HT/KT**par.prodalphaK[0])**(1/(1-par.prodalphaK[0]))  
            NNT[:] = (HNT/KNT**par.prodalphaK[1])**(1/(1-par.prodalphaK[1]))  
        else:
            NT[:] = (((HT/(par.TFP * par.TFP_s[0]))**exp_NK -par.prodalphaK[0] * KT**exp_NK)/par.prodbetaN[0])**(1/exp_NK) 
            NNT[:] = (((HNT/(par.TFP * par.TFP_s[1]))**exp_NK -par.prodalphaK[1] * KNT**exp_NK)/par.prodbetaN[1])**(1/exp_NK) 
    
        # Calculate marginalproduct of labor rWT
        rWT[:] = rHT * par.prodbetaN[0] * (par.TFP * par.TFP_s[0]) ** exp_NK * (HT / NT) ** (
                1 / par.sigma_NK)   
        rWNT[:] = rHNT * par.prodbetaN[1] * (par.TFP * par.TFP_s[1]) ** exp_NK * (HNT / NNT) ** (
                1 / par.sigma_NK)     
    else:
        NT[:] = HT / (par.TFP * par.TFP_s[0])
        NNT[:] = HNT / (par.TFP * par.TFP_s[1])
        rWT[:] = par.TFP_s[0]*par.TFP*rHT
        rWNT[:] = par.TFP_s[1]*par.TFP*rHNT

      
    # Wages from labor adjustment costs 
    adj_term_T = get_labor_adj_cost(r, ss.r, NT, ss.NT, par.phi_N)
    adj_term_NT = get_labor_adj_cost(r, ss.r, NNT, ss.NNT, par.phi_N)

    WT[:] = wT * P 
    WNT[:] = wNT * P 

    Lfoc_T[:] = WT*(1 - subP) - P*(rWT - adj_term_T)
    Lfoc_NT[:] = WNT*(1 - subP) - P*(rWNT - adj_term_NT)

    # Value added part of CES tree
    nOT = PH * ZT - PXT * XT - PH * par.FixedCost[0]
    nONT = PNT * ZNT - PXNT * XNT - PNT * par.FixedCost[1]
    N[:] = NT + NNT

    # GDP deflator 
    for t in range(par.T):
        nOT_denom = ss.PH * ZT[t] - ss.PXT * XT[t] - ss.PH * par.FixedCost[0]
        nONT_denom = ss.PNT * ZNT[t] - ss.PXNT * XNT[t] - ss.PNT * par.FixedCost[1]
        PGDP[t] = (nOT[t] + nONT[t]) / (nOT_denom + nONT_denom) * ss.PGDP
        GDP[t] = (nOT[t] + nONT[t]) / PGDP[t]

    YT[:] = nOT / PGDP
    YNT[:] = nONT / PGDP
    Y[:] = YT + YNT
    Z[:] = (PH * ZT + PNT * ZNT) / PGDP

    #####################
    #        Wages      #
    #####################
    # Average real wage income     
    wnT[:] = (WT * NT / P) / par.sT
    wnNT[:] = (WNT * NNT / P) / (1 - par.sT)

    # Aggregate real wage income
    labor_comp[:] = wT*NT + wNT*NNT 
    w[:] = labor_comp/N


@nb.njit
def asset_returns(par, ini, ss, DivT, DivNT, Div, YT, YNT, VAT, wnT,
                wnNT, subP, subI, P, PGDP, r, pD, ra, rai, AI_lag, AI, d, pi, IT, INT, 
                adjI_T, adjI_NT, mT, mNT, mT_res, mNT_res,
                KT, KNT, rKT, rKNT, Div_k_T, Div_k_NT):

    DivT[:] = (1 - VAT) * YT * PGDP / P - par.sT * wnT * (1 - subP) - IT * (1 - subI) 
    DivNT[:] = (1 - VAT) * YNT * PGDP / P  - (1 - par.sT) * wnNT * (1 - subP) - INT * (1 - subI) 
    Div[:] = DivT + DivNT
    Div_p = lead(Div, ss.Div)

    # Solve pD = (pD_p+Div)/(1+r) for pD
    sol_backwards_lin(pD, ss.pD, Div_p / (1 + r), 1 / (1 + r), par.T)

    if par.HH_type == 'HA-2A' or par.use1A_illiquid:
        tot_return = ((pD[0] + Div[0] + (1 + ss.i) / (1 + pi[0]) * ss.B) / ss.Atot - 1)
        rai[0] = par.xi + (tot_return*ss.Atot - ss.r * ss.A)/ss.AI 
        rai[1:] = par.xi + r[:-1]
        ra[:] = lag(ss.r, r)

        if par.use1A_illiquid:
            for t in range(par.T):
                AI_lag[t] = AI[t-1] if t > 0 else ss.AI 
                d[t] = ss.rai/(1+ss.rai)*(1+rai[t])*AI_lag[t] + par.chi*((1+rai[t])*AI_lag[t] -(1+ss.rai)*par.AI_target)
                AI[t] = (1+rai[t]) * AI_lag[t] - d[t]
            
    else:
        # Return at time 0 of the MIT shock is calculated from initial asset position
        ra[0] = (pD[0] + Div[0] + (1 + ss.i) / (1 + pi[0]) * ss.B) / ss.A - 1
        # For the remaining periods arbitrage-conditions hold so 
        ra[1:] = r[:-1]

        AI_lag[:] = 0.
        d[:] = 0. 

    # m residual   
    mT_res[:] = DivT - mT
    mNT_res[:] = DivNT - mNT



@nb.njit
def public_finances(par, ini, ss, VAT, YT, YNT, P, PGDP, subP, subI, 
                    WT, NT, WNT, NNT, UniformT_exo, UniformT, 
                    G_trans, G_exo, G_T, G_NT, G, PH, PNT, B, i, 
                    pi, tau,  G_budget, IT, INT):


    TR = VAT * (YT + YNT) * PGDP / P - subP * (WT * NT + WNT * NNT) / P - subI*(IT + INT)
    UniformT[:] = UniformT_exo
    G_trans[:] = G_exo
    G_T[:] = par.sGT_ss * ss.G + par.sGT * G_trans
    G_NT[:] = (1 - par.sGT_ss) * ss.G + (1 - par.sGT) * G_trans
    G[:] = (PH * G_T + PNT * G_NT) / P

    B_lag = lag(ss.B, B)
    i_lag = lag(ss.i, i)

    if par.debt_rule:  # Auclert et al debt rule

        tau[:] = ss.tau + par.epsB*(B_lag-ss.B)/ss.Y

    else:  # Complete debt financed policy for the first tauB periods

        if not isclose(ss.B,0.0):

            discont_rule(tau, par.T, par.tauB, par.deltaB, ss.tau, par.epsB, B, ss.B)

        else:

            tau[:] = ss.tau

    G_budget[:] = B + tau*(WT*NT+WNT*NNT)/P - (G + UniformT - TR + (1 + i_lag) / (1 + pi) * B_lag)

@nb.njit
def discont_rule(tau, T, tauB, deltaB, tau_ss, epsB, B, B_ss, gross=True):
    if gross:
        con = 0. 
    else:
        con = 1.

    for t in range(T):

        if t < tauB:
        
            tau[t] = tau_ss
        
        elif t <= tauB + deltaB:

            ttilde = (t - tauB) / deltaB
            tautilde = (con+tau_ss) * (B[t - 1] / B_ss) ** epsB
            omega_x = 3 * ttilde ** 2 - 2 * ttilde ** 3
            tau[t] = (1 - omega_x) * (con+tau_ss) + omega_x * tautilde - con

        else:

            tau[t] = (con+tau_ss) * (B[t] / B_ss) ** epsB - con

    

@nb.njit
def simple_HHs(par, ini, ss, wnNT, wnT, UniformT, tau, Income, C, r, 
                eps_beta, CHtM, CR, iF_s, piF, Q, C_s, 
                ra, A, Atot, UC_T, UC_NT):

    Income[:] = (1-tau)* (wnNT * (1 - par.sT) + wnT * par.sT) + UniformT
    
    if par.HH_type not in ['HA', 'HA-2A']:
        if par.use_RA_jac:
            r_lag = lag(ss.r, r)
            reval = (ra - r_lag)*ss.A 
            C[:] = ss.C
            C[:] += par.M_Y @ ((Income - ss.Income) + reval)
            C[:] += par.M_R @ (r - ss.r)
            C[:] += par.M_beta @ (eps_beta - ss.eps_beta) * par.beta_mean
        else:
            if par.HH_type == 'RA-CM' or par.HH_type == 'RA-IM':
                sHtM = 0. 
            else:
                sHtM = par.MPC_est[0]
            for tt in range(par.T):
                t = par.T - 1 - tt
                if par.HH_type == 'RA-CM' or par.HH_type == 'TA-CM':
                    if par.HH_type == 'TA-CM':
                        CHtM[t] = Income[t]
                    else:
                        CHtM[t] = 0.0
                    if t == par.T - 1:
                        CR[t] = ss.CR
                    else:
                        rF_s = iF_s[t] - piF[t]
                        CR[t] = CR[t + 1] * (Q[t + 1] / Q[t] * (1 + rF_s) / (1 + ss.r)) ** (-1 / par.CRRA) * C_s[t] / \
                                C_s[t + 1]
                    C[t] = CR[t] * (1 - sHtM) + sHtM * CHtM[t]
                elif par.HH_type == 'RA-IM' or par.HH_type == 'TA-IM':
                    if par.HH_type == 'TA-IM':
                        CHtM[t] = Income[t]
                    else:
                        CHtM[t] = 0.0
                    if t == par.T - 1:
                        CR[t] = ss.CR
                    else:
                        CR[t] = ((1 + ra[t + 1]) * eps_beta[t] * par.beta_mean * CR[t + 1] ** (-par.CRRA)) ** (
                                    -1 / par.CRRA)
                C[t] = CR[t] * (1 - sHtM) + sHtM * CHtM[t]

        UC_T[:] = C[:] ** (-par.CRRA)
        UC_NT[:] = C[:] ** (-par.CRRA)
        Get_HH_A(A, C, Income, ra, par.T, ss)
        Atot[:] = A


@nb.njit
def mon_pol(par, ini, ss, PNT, PH, DomP, ppi, pi, piNT, E, ND, Z, i, di, Taylor):

    i_lag = lag(ss.i, i)
    pi_p = lead(pi, ss.pi)

    # Calculate domestic price index using Paasche price index (sum of NT and H)
    DomP[:] = (PNT * ss.CNT + PH * ss.CH) / (ss.PNT * ss.CNT + ss.PH * ss.CH)
    ppi[:], ppi_p = Inf(DomP, ss.DomP)

    if par.floating:
        if par.MonPol == 'Taylor':
            if par.TaylorType == 'CPI':
                Taylor_pi = pi
            elif par.TaylorType == 'NT':
                Taylor_pi = lead(piNT, ss.piNT)
            elif par.TaylorType == 'ppi':
                Taylor_pi = ppi
            elif par.TaylorType == 'DomP':
                Taylor_pi = np.log(DomP) - np.log(ss.DomP)
            elif par.TaylorType == 'FoF':
                Taylor_pi = ppi_p + 0.3 / par.phi * np.log(E / lag(ss.E, E))
            elif par.TaylorType == 'Y':
                Taylor_pi = ppi_p + 0.25 / par.phi * np.log(Z / ss.Z)
            Taylor[:] = i - ((ss.i + par.phi * Taylor_pi) * (1 - par.phi_back) + par.phi_back * i_lag + di)
        else:
            Taylor[:] = 1 + ss.r - (1 + i + di) / (1 + pi_p)
    else:
        Taylor[:] = E - ND


@nb.njit
def pricing(par, ini, ss, PF_s, P, pi, piF_s, piF, PF, PH_s, piH_s,
            rF_s, iF_s, E, VAT, Q, r, i, NFA, PH, PT, PNT, 
            piNT, ToT, CH_s, C_s, piH, FD, subP, subI, pi_ann, piH_ann, 
            r_ann, rF_s_ann, piF_s_ann, iF_s_ann, r_NFA):

    # update parameters
    par.theta_w[:] = par.epsilon / par.NKWslope # warning: this is not thread safe
    par.theta[:] = par.epsilon / par.NKslope # warning: this is not thread safe

    # Fiscal devaluation - assumes that taxes are zero in steady state
    if par.FD_shock:
        subP[:] = par.subP_weight * FD
        subI[:] = par.subI_weight * FD
        VAT[:] = par.VAT_weight * FD
    else:
        subP[:] = 0.
        subI[:] = 0.
        VAT[:] = 0. 

    # Prices and inflation 
    piF_s_p = lead(piF_s, ss.piF_s)
    P_from_inf(PF_s, piF_s, par.T, ss.PF_s)
    P_from_inf(P, pi, par.T, ss.P)  # P
    P_from_inf(PF, piF, par.T, ss.PF)  # PF
    P_from_inf(PH_s, piH_s, par.T, ss.PH_s)  # PH_s
    pi_p = lead(pi, ss.pi)
    rF_s[:] = (1 + iF_s) / (1 + piF_s_p) - 1
    rF_s_ann[:] = convert_q2a(rF_s)
    piF_s_ann[:] = convert_q2a(piF_s)
    iF_s_ann[:] = convert_q2a(iF_s)

    # Exchange rate 
    E[:] = PF / PF_s * (1 - VAT)
    Q[:] = E * PF_s / P
    Q_terminal = ss.Q 
    Q_p = lead(Q, Q_terminal)

    # UIP  
    discont_rule(r_NFA, par.T, par.tauNFA, par.deltaNFA, ss.r_NFA, par.epsNFA, 1+NFA, 1., gross=False)
    NFA_closure = np.exp(-par.r_debt_elasticity * (NFA / ss.GDP - ss.NFA / ss.GDP)) - r_NFA
    if par.HH_type in ['HA', 'HA-2A']:
        port_adj_term = par.UIP_dev * (NFA/ss.Atot - lag(ss.NFA, NFA)/ss.Atot) - 1/(1+ss.r) * par.UIP_dev * (lead(NFA, ss.NFA)/ss.Atot - NFA/ss.Atot)
        r[:] = (1 + rF_s) * (Q_p / Q) * NFA_closure + port_adj_term - 1
    else:
        port_adj_term = par.UIP_dev * (NFA/ss.Atot - lag(ss.NFA, NFA)/ss.Atot) - 1/(1+ss.r) * par.UIP_dev * (lead(NFA, ss.NFA)/ss.Atot - NFA/ss.Atot)
        r[:] = (1 + rF_s) * (Q_p / Q) * NFA_closure + port_adj_term - 1
    
    # Fisher
    i[:] = (1 + r) * (1 + pi_p) - 1

    # LCP/PCP 
    PH[:] = PH_s * E / (1 - VAT)
    piH[:], piH_p = Inf(PH, ss.PH)

    PT[:] = Price_index(PH, PF, par.eta, par.alphaF)
    if isclose(par.etaT, 1.0):
        PNT[:] = (P / PT ** par.alphaT) ** (1 / ((1 - par.alphaT)))
    else:
        PNT[:] = ((P ** (1 - par.etaT) - par.alphaT * PT ** (1 - par.etaT)) / (1 - par.alphaT)) ** (1 / (1 - par.etaT))
    piNT[:], piNT_p = Inf(PNT, ss.PNT)

    ToT[:] = PF / PH_s / E * (1 - VAT)

    # Exports 
    CH_s[:] = Armington_dyn(PH_s, PF_s, C_s, par.gamma, par.alpha_s, par.phi_X, par.T, ss.PH_s, ss.PF_s)

    # Annual inflation and real rate    
    annual_inflation(pi_ann, P, par.T, ss.P)
    pi_ann_p = lead(pi_ann,0.)
    r_ann[:] = (1+i)**4 / (1+pi_ann_p) - 1
    annual_inflation(piH_ann, PH, par.T, ss.PH)


@nb.njit(cache=False)
def block_post(par, ini, ss, P, CH_s, PT, PF, PNT, PH, 
                ZT, ZNT, GDP_NT, GDP_T, pi,  
                XT, XNT, XNT2T, XT2NT, XNT2NT, 
                XT2T, XM2T, PXT, PXNT,
                NT, NNT, WNT, WT, wNT, wT, piWT, piWNT,
                pD, ra, G, G_T, G_NT, B,
                UC_T, UC_NT, C, C_T, C_NT, A, CT, CNT, CF,
                CH, NFA, NX, Exports, Imports, Walras,
                IT, INT, I, adjI_T, adjI_NT, 
                NFA_target, NKWPCNT, NKWPCT, 
                goods_mkt_NT, goods_mkt_T, 
                C_T_hh, C_NT_hh, C_hh, A_hh, AI_hh, AI, Atot, UC_T_hh, UC_NT_hh, tau):

    if par.HH_type == 'HA':
        C_T[:] = C_T_hh
        C_NT[:] = C_NT_hh
        C[:] = C_hh
        A[:] = A_hh
        if par.use1A_illiquid:
            Atot[:] = AI + A
        else:
            Atot[:] =  A

        UC_T[:] = UC_T_hh
        UC_NT[:] = UC_NT_hh

    elif par.HH_type == 'HA-2A':
        C_T[:] = C_T_hh
        C_NT[:] = C_NT_hh
        AI[:] = AI_hh
        C[:] = C_hh  
        A[:] = A_hh
        Atot[:] = AI + A

        UC_T[:] = UC_T_hh
        UC_NT[:] = UC_NT_hh   
    else:
        C_T[:] = C
        C_NT[:] = C

    # Get consumption at lower levels of CES tree
    CT[:], CNT[:] = CES_demand_T(C, PT, PNT, P, par.etaT, par.alphaT)
    CF[:], CH[:] = CES_demand(CT, PH, PF, PT, par.eta, par.alphaF)

    ######################
    #       Walras       #
    ######################        
    NFA_target[:] = NFA - (Atot - pD - B)
    NFA_lag = lag(ss.NFA, NFA)
    NX[:] = (PH * ZT + PNT * ZNT - PXT * XT - PXNT * XNT) / P - C - G - IT  - INT  - PH / P * par.FixedCost[0] - PNT / P * \
            par.FixedCost[1]
    Walras[:] = NFA - (NX + (1 + ra) * NFA_lag)
    Exports[:] = CH_s
    Imports[:] = CF + XM2T

    ####################################
    # Goods market clearing and NKWPCs #
    ####################################

    GDP_T[:] = ((ZT - par.FixedCost[0]) - XT2NT - XT2T)
    GDP_NT[:] = ((ZNT - par.FixedCost[1]) - XNT2T - XNT2NT)
    goods_mkt_T[:] = GDP_T * PH - (PH * CH + PH * CH_s + PH * G_T + P * IT)
    goods_mkt_NT[:] = GDP_NT * PNT - (PNT * CNT + PNT * G_NT + P * INT)

    PC_beta = 1 / (1 + ss.r)
    piWT[:], piWT_p = Inf(WT, ss.WT)
    piWNT[:], piWNT_p = Inf(WNT, ss.WNT)
    pi_l = lag(ss.pi, pi)

    v_prime = par.psi[0] * (NT / par.sT) ** par.inv_frisch
    LHS_WPC_T = ((1+piWT)/(1+pi_l)**par.piW_index - 1) * (1+piWT)/(1+pi_l)**par.piW_index 
    RHS_WPC_T = par.NKWslope[0] * NT * (v_prime - (1-tau)*wT * 1/par.muw[0] * UC_T) + PC_beta * ((1+piWT_p)/(1+pi)**par.piW_index - 1) * (1+piWT_p)/(1+pi)**par.piW_index
    
    NKWPCT[:] = LHS_WPC_T - RHS_WPC_T

    v_prime = par.psi[1] * (NNT / (1 - par.sT)) ** par.inv_frisch
    LHS_WPC_NT = ((1+piWNT)/(1+pi_l)**par.piW_index - 1) * (1+piWNT)/(1+pi_l)**par.piW_index
    RHS_WPC_NT = par.NKWslope[1] * NNT * (
                v_prime - (1-tau)*wNT * 1/par.muw[1] * UC_NT) + PC_beta * ((1+piWNT_p)/(1+pi)**par.piW_index - 1) * (1+piWNT_p)/(1+pi)**par.piW_index
    

    NKWPCNT[:] = LHS_WPC_NT - RHS_WPC_NT

