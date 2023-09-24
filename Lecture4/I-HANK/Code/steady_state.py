import time
import numpy as np
from scipy import optimize
from scipy.optimize import leastsq

from consav.misc import elapsed

from broyden_solver_autojac import broyden_solver_autojac
from household_problem import compute_RA_jacs
import blocks

class StopExecution(Exception): pass # used in calibration

def calibrate_ss(model, do_print=False, print_timeings=False):
    """ find the steady state """

    t0 = time.time()

    par = model.par
    ss = model.ss

    # Prices normalized to 1
    P_list = ['E', 'PNT', 'PH', 'PF_s', 'PT', 'PXT', 'PXNT', 'PH']
    for i in P_list:
        setattr(ss, i, 1.0)

    # Inflation rates are zero in the steady state
    pi_list = ['piNT', 'pi', 'piF', 'piH', 'ppi', 'piF_s', 'piH_s', 'piWT', 'piWNT']
    for i in pi_list:
        setattr(ss, i, 0.0)

    # Use model relations from here on
    ss.PH_s, ss.PF = ss.PH / ss.E, ss.PF_s * ss.E
    ss.PTP = ss.PT / ss.P
    ss.P = blocks.Price_index(ss.PNT, ss.PT, par.etaT, par.alphaT)
    ss.Q = ss.E / ss.P

    # Normalize aggregate GDP and labor to 1
    ss.GDP = 1.0
    ss.N = 1.0

    # Nominal interest rates are all equal to foreign, exogenous rate
    ss.i = ss.iF_s = ss.iF_s = ss.rF_s = par.iF_s_exo

    # Real rates are equal to nominal rates due to zero inflation
    ss.r = ss.ra = ss.i

    # Debt elasticity 
    ss.r_NFA = 0. 

    # Annual rates 
    ss.pi_ann = 0. 
    ss.piH_ann = 0.
    ss.r_ann = (1+ss.i)**4 / (1+ss.pi_ann) - 1
    ss.rF_s_ann = (1+ss.rF_s)**4 - 1
    ss.iF_s_ann = ss.rF_s_ann
    ss.piF_s_ann = 0.

    # shocks
    ss.C_s = 1.0
    ss.Y_s = 1.0 
    ss.eps_beta = 1.0
    ss.di = 0.0
    par.TFP = 1.0

    # Public sector   
    ss.G_exo = ss.G_trans = 0.0
    ss.UniformT = ss.UniformT_exo = 0.0
    ss.VAT = 0.0
    ss.subP = 0.0
    ss.subI = 0.0 
    ss.FD = 0.0

    # Preliminaries 
    ss.mcT = 1 / par.mu[0]
    ss.mcNT = 1 / par.mu[1]
    ss.PH_s = ss.PH / ss.E * (1 - ss.VAT)
    ss.PF = ss.E * ss.PF_s / (1 - ss.VAT)
    ss.ToT = ss.PF / ss.PH_s / ss.E * (1 - ss.VAT)
    ss.ND = ss.E

    par.epsilon[:] = par.mu/(par.mu-1)
    par.epsilon_w[:] = par.muw/(par.muw-1)
    par.theta[:] = par.epsilon / par.NKslope
    par.theta_w[:] = par.epsilon_w / par.NKWslope

    # define steady state residuals
    def ss_res(x):

        par.FixedCost[0], par.FixedCost[1], ss.ZNT, ss.ZT, par.X_share_T[1] = x[0], x[1], x[2], x[3], x[4]
        
        # Last element of X_share_T we get residually since sum(X_share_T) = 1 
        par.X_share_T[2] = 1 - par.X_share_T[0] - par.X_share_T[1]

        ss.NT = ss.N * par.sT
        ss.NNT = ss.N * (1 - par.sT)

        # solve supply side
        if par.X_expshare[0] > 0.0 and par.X_expshare[1] > 0.0:

            # solve first nest (X vs H)
            solve_production_XH(ss, par)

            if par.use_capital:
                
                solve_production_KN(ss, par)

            else:

                solve_production_N(ss, par)
                ss.KT, ss.KNT = 0., 0.
                ss.IT, ss.INT = 0., 0.
                ss.rKT, ss.rKNT = 0., 0.
                ss.I = 0.
                ss.ITp, ss.INTp = 0., 0.
                ss.KTp, ss.KNTp = 0., 0.
                ss.adjI_T = 0.
                ss.adjI_NT = 0.

        else:
            ss.XT = 0. 
            ss.XNT = 0. 
            ss.HT = ss.ZT 
            ss.HNT = ss.ZNT 
            ss.rHT = (1 - ss.VAT) * ss.mcT * ss.PH  / ss.P
            ss.rHNT = (1 - ss.VAT) * ss.mcNT * ss.PNT / ss.P

            # If we have capital and labor 
            if par.use_capital:

                solve_production_KN(ss, par)

            # If we do not have intermediate goods (X_expshare=0) or capital but only labor in production 
            else:

                solve_production_N(ss, par)

        # Get demand for intermediate goods from different sectors using CES formulas
        ss.XT2NT, ss.XNT2T, ss.XNT2NT, ss.XT2T, ss.XM2T = blocks.get_intermediates(ss.XT, ss.XNT, ss.PXT, ss.PXNT,
                                                                                    ss.PH, ss.PNT, ss.PF, par.X_share_T,
                                                                                    par.X_share_NT, par.etaX)

        # Nominal value added is output net intermediate goods and fixed cost
        nOT = ss.PH * ss.ZT - ss.PXT * ss.XT - ss.PH * par.FixedCost[0]
        nONT = ss.PNT * ss.ZNT - ss.PXNT * ss.XNT - ss.PNT * par.FixedCost[1]
        ss.PGDP = (nOT + nONT) / ss.GDP  # GDP deflator is nominal GDP over real GDP

        # Real value added is nominal value added deflated with the GDP deflator
        ss.YT = nOT / ss.PGDP
        ss.YNT = nONT / ss.PGDP

        # Other aggregates and accounting
        ss.Y = ss.YT + ss.YNT
        ss.Z = (ss.PH * ss.ZT + ss.PNT * ss.ZNT) / ss.PGDP
        ss.PO_T = nOT / ss.GDP_T
        ss.PO_NT = nONT / ss.GDP_NT
        ss.OT = ss.YT * ss.PGDP / ss.PO_T
        ss.ONT = ss.YNT * ss.PGDP / ss.PO_NT

        # Wages
        ss.wnT = (ss.WT * ss.NT / ss.P) / par.sT
        ss.wnNT = ((ss.WNT * ss.NNT / ss.P) / (1 - par.sT))

        ss.wT, ss.wNT = ss.wnT * par.sT / ss.NT, ss.wnNT * (1 - par.sT) / ss.NNT
        ss.labor_comp = ss.wT*ss.NT + ss.wNT*ss.NNT 
        ss.w = ss.labor_comp / ss.N 
        
        # Dividends and asset pricing
        ss.DivT = ((1 - ss.VAT) * nOT - ss.WT * ss.NT * (1 - ss.subP)) / ss.P - ss.IT
        ss.DivNT = ((1 - ss.VAT) * nONT - ss.WNT * ss.NNT * (1 - ss.subP)) / ss.P - ss.INT
        
        if par.use_capital:
            ss.Div_k_T = ss.rKT*ss.KT - ss.IT 
            ss.Div_k_NT = ss.rKNT*ss.KNT - ss.INT 
            ss.RP_T  = 1.
            ss.RP_NT  = 1.

        ss.Div = ss.DivT + ss.DivNT
        ss.pD = ss.Div / ss.r

        # Public sector (bonds, taxes, spending)
        ss.G = par.G_GDP_ratio * ss.GDP * ss.PGDP / ss.P
        ss.B = par.B_GDP_ratio * ss.GDP * ss.PGDP / ss.P
        Tax_revenue = ss.VAT * (nOT + nONT) / ss.P - ss.subP * (ss.WT * ss.NT + ss.WNT * ss.NNT) / ss.P
        ss.tau = (ss.r * ss.B + ss.G - Tax_revenue)/((ss.WT * ss.NT + ss.WNT * ss.NNT)/ss.P)

        # Households
        ss.Income = (1-ss.tau)*(ss.wnNT * (1 - par.sT) + ss.wnT * par.sT) + ss.UniformT + ss.UniformT_exo
        A = ss.pD + ss.B  # or A = par.W2INC_target * ss.Income
        ss.C = ss.r * A + ss.Income

        # Aggregate trade flows

        # Exports
        nExports = par.X2Y_target * ss.GDP
        ss.Exports = nExports * ss.PGDP / ss.P

        # Using that imports = exports in steady state (NFA = 0 by assumption)
        ss.Imports = ss.Exports
        HH_imports = ss.Imports - ss.XM2T * ss.PF / ss.P
    
        # Back out alpha that is consistent with implied level of HH imports
        par.alphaF = HH_imports / (par.alphaT * ss.C)

        # Get consumption of various goods from CES
        ss.CT, ss.CNT = blocks.CES_demand(ss.C, ss.PNT, ss.PT, ss.P, par.etaT, par.alphaT)
        ss.CF, ss.CH = blocks.CES_demand(ss.CT, ss.PH, ss.PF, ss.PT, par.eta, par.alphaF)
        
        # Back out alphaF such that foreign demand = exports
        par.alpha_s = ss.Exports * ss.PGDP / ss.P
        ss.CH_s = blocks.Armington(ss.PH_s, ss.PF_s, ss.C_s, par.gamma, par.alpha_s)

        # Goods market clearing
        ss.GDP_T = ss.ZT - ss.XT2NT - ss.XT2T - par.FixedCost[0]
        ss.GDP_NT = ss.ZNT - ss.XNT2T - ss.XNT2NT - par.FixedCost[1]
        ss.goods_mkt_T = ss.GDP_T - (ss.CH + ss.CH_s + ss.G * par.sGT_ss + ss.IT)
        ss.goods_mkt_NT = ss.GDP_NT - (ss.CNT + ss.G * (1 - par.sGT_ss) + ss.INT)

        # Calibration targets
        Asset_target = (ss.pD + ss.B) / ss.Income - par.W2INC_target  # wealth-income ratio
        profit_target = ss.DivT / ss.Div - ss.GDP_T / ss.GDP  # Dividend share in sectors should be same as GDP as share

        if par.X_expshare[0] > 0.0 and par.X_expshare[1] > 0.0:
            Import_Target = HH_imports / ss.Imports - par.HH_importshare
        else:
            Import_Target = par.X_share_T[2] - 0.0

        return np.array([Asset_target, ss.goods_mkt_NT, ss.goods_mkt_T, profit_target, Import_Target])

    # Calibrate steady state
    if do_print: print('Solving calibration:')

    results = broyden_solver_autojac(ss_res,x0=par.x0,maxiter=150,tol=1E-10)
    ss_res(results[0])  # evaluation at solution

    par.x0 = results[0]

    assert np.isclose(sum(par.X_share_T), 1.0)
    assert np.min(par.X_share_T) >= 0
    assert par.alpha_s >= 0
    assert par.alphaF >= 0
    assert par.alphaF <= 1

    if print_timeings: print(f'Step 1 in calibration solved in {elapsed(t0)}')

    #########################
    #   Household problem   #
    #########################

    t1 = time.time()

    if par.HH_type in ['HA', 'HA-2A']:
    
        if par.HH_type == 'HA':
    
            if par.use1A_illiquid:
    
                # Return to liquid asset 
                ss.rai = ss.ra + par.xi
                par.rai_ss_target = ss.rai

                def HH_res(x):

                    par.beta_delta = 0. # We only have one beta with HA illiquid 
                    par.beta_mean, par.AI_target  = x

                    ss.AI = par.AI_target
                    ss.AI_lag = ss.AI 
                    ss.d = ss.rai/(1+ss.rai)*(1+ss.rai)*ss.AI+ par.chi*((1+ss.rai)*ss.AI -(1+ss.rai)*par.AI_target)                    

                    model.solve_hh_ss(do_print=False)
                    model.simulate_hh_ss(do_print=False)

                    ss.A = np.sum(ss.D * ss.a)
                    ss.C = np.sum(ss.D * ss.c) 

                    MPC_ann, MPC_quarterly = model.ann_MPCs()
                    
                    MPC_y1_res = MPC_ann[0] - par.MPC_est[0]
                    MPC_y2_res = MPC_ann[1] - par.MPC_est[1]
                    
                    A_res = ss.A + ss.AI - (ss.pD + ss.B)

                    if np.abs(A_res) < 1E-11 and np.abs(MPC_y1_res) < 1E-11:
                        raise StopIteration

                    return np.array([A_res, MPC_y1_res])

                if do_print: print('Solving HA calibration:')

                from scipy import optimize
                
                try:
                    
                    results_HH = optimize.root(HH_res, np.array([0.98146609, 5.02143341]),  method='lm')

                    par.x0_het = results_HH.x
                    par.beta_mean = results_HH.x[0]
                    par.AI_target = results_HH.x[1]
                    HH_res(results_HH.x)

                except StopIteration:

                    pass
                
            else:

                ss.AI_lag = 0. 
                ss.d = 0.

                def HH_res(x):

                    par.beta_mean, par.beta_delta = x[0], x[1]

                    maxbeta = max(par.beta_grid)
                    if maxbeta >= 0.998:
                        print('Max beta reached!')
                        maxbeta = 0.9998 * 1 / (1 + ss.r) - 0.0005
                        par.beta_delta = maxbeta - par.beta_mean

                    model.solve_hh_ss(do_print=False)
                    model.simulate_hh_ss(do_print=False)

                    ss.A = np.sum(ss.D * ss.a)
                    ss.C = np.sum(ss.D * ss.c)

                    Target_first_year_MPC = True

                    MPC_ann, MPC_quarterly = model.ann_MPCs()
                    if Target_first_year_MPC:
                        MPC_res = MPC_ann[0] - par.MPC_est[0]
                    else:
                        MPC_res = MPC_quarterly[0] - 0.19961746
                    
                    A_res = (ss.A - ss.pD - ss.B)
                    if np.abs(A_res) < 1E-11 and np.abs(MPC_res) < 1E-11:
                        raise StopIteration

                    return np.array([A_res, MPC_res])

                if do_print: print('Solving HA calibration:')

                try:

                    results_HH = broyden_solver_autojac(HH_res, x0=par.x0_het,maxiter=30,tol=1E-11)
                    par.x0_het = results_HH[0]
                    HH_res(results_HH[0])
            
                except StopIteration:

                    pass

        elif par.HH_type == 'HA-2A':

            # Return to liquid asset 
            ss.rai = ss.ra + par.xi
            par.rai_ss_target = ss.rai

            def HH_res(x):

                par.beta_delta = 0. # We only have one beta with HA-2A 
                par.beta_mean, par.AI_target  = x

                model.solve_hh_ss(do_print=False)
                model.simulate_hh_ss(do_print=False)

                ss.A = np.sum(ss.D * ss.a)
                ss.AI = np.sum(ss.D * ss.ai)
                ss.C = np.sum(ss.D * ss.c)  

                MPC_ann, MPC_quarterly = model.ann_MPCs()
                
                MPC_y1_res = MPC_ann[0] - par.MPC_est[0]
                MPC_y2_res = MPC_ann[1] - par.MPC_est[1]

                A_res = ss.A + ss.AI - (ss.pD + ss.B)

                if np.abs(A_res) < 1E-11 and np.abs(MPC_y1_res) < 1E-11:
                    raise StopIteration

                return np.array([A_res, MPC_y1_res])

            if do_print: print('Solving HA-2A calibration:')

            try:

                results_HH = broyden_solver_autojac(HH_res,x0=par.x0_het,maxiter=30,tol=1E-11)
            
                par.x0_het = results_HH[0]
                par.beta_mean = results_HH[0][0]
                par.AI_target = results_HH[0][1]
                HH_res(results_HH[0])

            except StopIteration:

                pass

        ss.A = np.sum(ss.D * ss.a)
        ss.C = np.sum(ss.D * ss.c) 
        ss.C_T = np.sum(ss.D * (ss.c_T ))
        ss.C_NT = np.sum(ss.D * (ss.c_NT))
        ss.C_T_hh = np.sum(ss.D * ss.c_T)
        ss.C_NT_hh = np.sum(ss.D * ss.c_NT)
        ss.UC_T = ss.UC_T_hh
        ss.UC_NT = ss.UC_NT_hh

        ss.C_T = np.sum(ss.D * ss.c_T)
        ss.C_T = np.sum(ss.D * ss.c_T)

        ss.C_hh = ss.C  
        ss.A_hh = ss.A

        ss.INC_T_hh = np.sum(ss.D * ss.inc_T)
        ss.INC_NT_hh = np.sum(ss.D * ss.inc_NT)

        if par.HH_type == 'HA':

            if par.use1A_illiquid: 
                ss.Atot = ss.A + ss.AI 
            else:
                ss.Atot = ss.A 

        elif par.HH_type == 'HA-2A':

            ss.Atot = ss.A + ss.AI 

    elif par.HH_type in par.No_HA_list:

        ss.A = ss.pD + ss.B
        ss.Atot = ss.A
        par.beta_mean = 1 / (1 + ss.r)
        ss.UC_T, ss.UC_NT = ss.C ** (-par.CRRA), ss.C ** (-par.CRRA)
        ss.C_T = ss.C_NT = ss.C

        if par.HH_type == 'TA-CM' or par.HH_type == 'TA-IM':
            ss.CHtM = ss.Income
            ss.CR = ss.r * ss.A / (1 - par.MPC_est[0]) + ss.Income            
        else:
            ss.CR = ss.C
            ss.CHtM = 0.0

        compute_RA_jacs(ss, par)

    else:

        raise ValueError('Incorrect HH type chosen!')

    if print_timeings:

        print(f'Step 2 in calibration solved in {elapsed(t1)}')

    # Back out disutility of work psi such that NKPWC holds in steady state
    par.psi[0] = (1-ss.tau)*ss.wT * 1/par.muw[0] * ss.UC_T / (ss.NT / par.sT) ** par.inv_frisch
    par.psi[1] = (1-ss.tau)*ss.wNT * 1/par.muw[1] * ss.UC_NT / (ss.NNT / (1 - par.sT)) ** par.inv_frisch

    # misc
    ss.NFA = ss.Atot - ss.pD - ss.B
    ss.goods_mkt = ss.Z * ss.PGDP - ss.XT2NT - ss.XNT2T - ss.XNT2NT - ss.XT2T - par.FixedCost[0] - par.FixedCost[1] - (
                ss.CH_s + ss.CNT + ss.CH + ss.G + ss.I)
    ss.goods_mkt_T = ss.ZT - ss.XT2NT - ss.XT2T - par.FixedCost[0] - (ss.CH + ss.CH_s + ss.G * par.sGT_ss + ss.IT)
    ss.goods_mkt_NT = ss.ZNT - ss.XNT2T - ss.XNT2NT - par.FixedCost[1] - (ss.CNT + ss.G * (1 - par.sGT_ss)  + ss.INT)

    assert np.isclose(ss.Z * ss.PGDP, ss.ZT * ss.PH + ss.ZNT * ss.PNT)
    assert np.isclose(ss.goods_mkt, 0., atol=1e-07)

    ss.NX = (ss.PH * ss.ZT + ss.PNT * ss.ZNT) / ss.P - (
                ss.C + ss.I + ss.G + ss.XT2NT + ss.XNT2T + ss.XNT2NT + ss.XT2T + ss.XM2T) - par.FixedCost[0] - par.FixedCost[1]  # + r_pay
    ss.Exports = ss.CH_s * ss.PH_s / ss.P
    ss.Imports = (ss.CF + ss.XM2T) * ss.PF / ss.P
    assert np.isclose(ss.Exports - ss.Imports, ss.NX, atol=1e-07)

    assert abs(ss.C - ((1-ss.tau)*(ss.WNT * ss.NNT + ss.WT * ss.NT) + ss.r * ss.Atot)) < 1e-07
    ss.Walras = ss.NX + ss.r * ss.NFA
    ss.DomP = (ss.PNT * ss.CNT + ss.PH * ss.CH) / (ss.PNT * ss.CNT + ss.PH * ss.CH)

    # mT and mNT if capital is used
    if par.use_capital:
        ss.mT = ss.DivT
        ss.mNT = ss.DivNT
        
    assert np.isclose(ss.rKT, ss.rKNT)

def solve_production_XH(ss, par):

    # Intermediate good FOCs
    exp = (par.sigma_XH - 1) / par.sigma_XH
    if par.X_expshare[0] > 0.0 and par.X_expshare[1] > 0.0:

        def Intermediate_goods_sys(x):

            par.alphaX[0], par.alphaX[1], ss.XT, ss.XNT = x
            exp = (par.sigma_XH - 1) / par.sigma_XH

            # Impose Cobb-douglas for identification of TFP 
            ss.HT = (ss.ZT/ (ss.XT ** par.alphaX[0] ))**(1/(1 - par.alphaX[0])) 
            ss.HNT = (ss.ZNT/ (ss.XNT ** par.alphaX[1]))**(1/(1 - par.alphaX[1])) 

            beta_T = (1 - par.alphaX[0]) * (ss.ZT / (ss.HT )) ** exp
            beta_NT = (1 - par.alphaX[1]) * (ss.ZNT / (ss.HNT)) ** exp

            ss.rHT = (1 - ss.VAT) * ss.mcT * beta_T * ss.PH * (ss.ZT / ss.HT) ** (
                    1 / par.sigma_XH)  / ss.P
            ss.rHNT = (1 - ss.VAT) * ss.mcNT * beta_NT * ss.PNT * (
                    ss.ZNT / ss.HNT) ** (1 / par.sigma_XH)  / ss.P

            # Nominal rental rates of labor 
            n_rHT = ss.rHT * ss.P
            n_rHNT = ss.rHNT * ss.P

            alpha_T = par.alphaX[0] * (ss.ZT / (ss.XT)) ** exp
            alpha_NT = par.alphaX[1] * (ss.ZNT / (ss.XNT)) ** exp

            res1 = ss.XT - alpha_T ** (par.sigma_XH) * (ss.PXT / ((1 - ss.VAT) * ss.PH * ss.mcT)) ** (
                -par.sigma_XH) * ss.ZT 
            res2 = ss.XNT - alpha_NT ** (par.sigma_XH) * (ss.PXNT / ((1 - ss.VAT) * ss.PNT * ss.mcNT)) ** (
                -par.sigma_XH) * ss.ZNT 

            res3 = ss.XT - par.X_expshare[0] * ss.HT * n_rHT  / (ss.PXT * (1 - par.X_expshare[0]))
            res4 = ss.XNT - par.X_expshare[1] * ss.HNT * n_rHNT  / (ss.PXNT * (1 - par.X_expshare[1]))

            return [res1, res2, res3, res4]

    solution = optimize.root(Intermediate_goods_sys, [0.5, 0.5, 1.0, 1.0], method='lm', options={'ftol': 1e-08})
    
    if not solution.success:
        raise ValueError('Could not solve intermediate goods system')
    
    residuals = Intermediate_goods_sys(solution.x)

    # check that solution is correct
    for k in residuals:
        assert abs(k) < 1e-06

    # get solution and update parameters
    par.alphaX[0], par.alphaX[1], ss.XT, ss.XNT = solution.x
    par.prodalphaX[0] = par.alphaX[0] * (ss.ZT / (ss.XT)) ** exp
    par.prodalphaX[1] = par.alphaX[1] * (ss.ZNT / (ss.XNT)) ** exp
    par.prodbeta[0] = (1 - par.alphaX[0]) * (ss.ZT / (ss.HT)) ** exp
    par.prodbeta[1] = (1 - par.alphaX[1]) * (ss.ZNT / (ss.HNT)) ** exp

def solve_production_KN(ss, par):

    # Intermediate good FOCs
    exp = (par.sigma_NK - 1) / par.sigma_NK
    ss.rKT  = ss.r + par.delta_K
    ss.rKNT = ss.r + par.delta_K
    def Intermediate_goods_sys(x):

        par.alphaK[0], par.alphaK[1],  ss.KT, ss.KNT, par.TFP_s[0], par.TFP_s[1] = x

        exp = (par.sigma_NK - 1) / par.sigma_NK
        NT_eff = ss.NT
        NNT_eff = ss.NNT

        beta_T = (1 - par.alphaK[0]) * (ss.HT / (NT_eff * par.TFP * par.TFP_s[0])) ** exp
        beta_NT = (1 - par.alphaK[1]) * (ss.HNT / (NNT_eff * par.TFP * par.TFP_s[1])) ** exp

        # capital shares
        alphaK_T = par.alphaK[0] * (ss.HT / (ss.KT * par.TFP * par.TFP_s[0])) ** exp
        alphaK_NT = par.alphaK[1] * (ss.HNT / (ss.KNT * par.TFP * par.TFP_s[1])) ** exp

        ss.rWT = ss.rHT * beta_T * (par.TFP * par.TFP_s[0]) ** exp * (ss.HT / NT_eff) ** (
                1 / par.sigma_NK) / (1 - ss.subP) / ss.P
        ss.rWNT = ss.rHNT * beta_NT * (par.TFP * par.TFP_s[1]) ** exp * (
                ss.HNT / NNT_eff) ** (1 / par.sigma_NK) / (1 - ss.subP) / ss.P

        # Nominal rental rates of labor 
        n_rWT = ss.rWT * ss.P
        n_rWNT = ss.rWNT * ss.P

        labor_comp_T = ss.NT * n_rWT * (1 - ss.subP)
        labor_comp_NT = ss.NNT * n_rWNT * (1 - ss.subP)
        capital_comp_T = ss.rKT * ss.KT
        capital_comp_NT = ss.rKNT * ss.KNT

        res1 = ss.KT - alphaK_T ** (par.sigma_NK) * (ss.rKT / ss.rHT) ** (
            -par.sigma_NK) * ss.HT * (par.TFP * par.TFP_s[0]) ** (par.sigma_NK - 1)
        res2 = ss.KNT - alphaK_NT ** (par.sigma_NK) * (ss.rKNT / ss.rHNT) ** (
            -par.sigma_NK) * ss.HNT * (par.TFP * par.TFP_s[1]) ** (par.sigma_NK - 1)

        res3 = ss.HT - (par.TFP * par.TFP_s[0] * NT_eff ** (1 - par.alphaK[0]) * ss.KT ** par.alphaK[0])
        res4 = ss.HNT - (par.TFP * par.TFP_s[1] * NNT_eff ** (1 - par.alphaK[1]) * ss.KNT ** par.alphaK[1])

        # Expenditure on capital
        res5 = par.K_expshare[0] - capital_comp_T/(capital_comp_T+labor_comp_T+ss.PXT*ss.XT)
        if np.isclose(par.K_expshare[1], 0.):
            res6 = par.alphaK[1] - 0.
        else:
            res6 = par.K_expshare[1] - capital_comp_NT / (capital_comp_NT + labor_comp_NT + ss.PXNT * ss.XNT)

        return [res1, res2, res3, res4, res5, res6]

    # alphaKT_guess, alphaKNT_guess = par.K_expshare
    alphaKT_guess, alphaKNT_guess = 0.3, 0.3

    solution = optimize.root(Intermediate_goods_sys, np.array([alphaKT_guess, alphaKNT_guess, 7.0, 10.0, 4.0, 1.5]), method='lm', options={'ftol': 1e-08})

    if not solution.success:
        raise ValueError('Could not solve intermediate goods system')
    
    residuals = Intermediate_goods_sys(solution.x)

    # check that solution is correct
    for k in residuals:
        assert abs(k) < 1e-06

    # get solution and update parameters
    par.alphaK[0], par.alphaK[1], ss.KT, ss.KNT, par.TFP_s[0], par.TFP_s[1] = solution.x
    par.prodalphaK[0] = par.alphaK[0] * (ss.HT / (ss.KT * par.TFP * par.TFP_s[0])) ** exp
    par.prodalphaK[1] = par.alphaK[1] * (ss.HNT / (ss.KNT * par.TFP * par.TFP_s[1])) ** exp if par.K_expshare[1]> 0  else 0.
    par.prodbetaN[0] = (1 - par.alphaK[0]) * (ss.HT / (ss.NT * par.TFP * par.TFP_s[0])) ** exp
    par.prodbetaN[1] = (1 - par.alphaK[1]) * (ss.HNT / (ss.NNT * par.TFP * par.TFP_s[1])) ** exp

    if not np.isclose(par.sigma_NK,
                      1.0):  # If we do not have Cobb-Douglass, check that CES is calibrated to match CD form in steady state
        assert np.isclose(ss.HT, par.TFP * par.TFP_s[0] * (
                par.prodbetaN[0] * ss.NT ** exp + par.prodalphaK[0] * ss.KT ** exp) ** (1 / exp))
        assert np.isclose(ss.HNT, par.TFP * par.TFP_s[1] * (
                par.prodbetaN[1] * ss.NNT ** exp + par.prodalphaK[1] * ss.KNT ** exp) ** (1 / exp))

    # Wages 
    ss.WT = ss.rWT * ss.P
    ss.WNT = ss.rWNT * ss.P

    # investment
    ss.IT = ss.KT * par.delta_K
    ss.INT = ss.KNT * par.delta_K
    ss.I = ss.IT + ss.INT
    ss.QT = 1.
    ss.QNT = 1.
    ss.ITp = ss.IT
    ss.INTp = ss.INT
    ss.KTp = ss.KT
    ss.KNTp = ss.KNT   
    ss.adjI_T = 0.
    ss.adjI_NT = 0.

def solve_production_N(ss, par):

    NT_eff = ss.NT
    NNT_eff = ss.NNT
    par.TFP_s[0] = ss.HT / (par.TFP * NT_eff)
    par.TFP_s[1] = ss.HNT / (par.TFP * NNT_eff)
    ss.WT = par.TFP_s[0]*par.TFP*ss.rHT
    ss.WNT = par.TFP_s[1]*par.TFP*ss.rHNT
    ss.rWT = ss.WT/ss.P
    ss.rWNT = ss.WNT/ss.P