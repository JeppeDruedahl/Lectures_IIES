import os
import time
import pickle
from copy import deepcopy
from types import SimpleNamespace
from IPython.display import display

import numpy as np
import pandas as pd
from scipy import ndimage,optimize

from EconModel import EconModelClass
from GEModelTools import GEModelClass

import steady_state
import household_problem

class StopIteration(Exception): pass # used in estimation

def create_model(name,par=None,do_print=True,ss=True,jac=True,
                 find_transition_path=True,skip_shocks=True,skip_hh=False):
    """ create model and solve for transition path to foreign shock """

    if do_print: print(name)
    if par is None: par = {}

    if not 'HH_type' in par:
        model = HANKModelClass(name,par=par)
        model.load_est_par(f'HANK_est_par_best')
        model.load_x0s(f'HANK_x0s')
    else:
        if par['HH_type'] == 'HA':
            model = HANKModelClass(name,par=par)
            model.load_est_par(f'HANK_est_par_best')
            model.load_x0s(f'HANK_x0s')
        elif par['HH_type'] in ['RA-IM','RA-CM']:
            model = HANKModelClass_RA(name,par=par)
            model.load_est_par(f'RANK_est_par_best')
        elif par['HH_type'] in ['TA-IM','TA-CM']:
            model = HANKModelClass_RA(name,par=par)
            model.load_est_par(f'HANK_est_par_best')
        elif par['HH_type'] in ['HA-2A']:
            model = HANKModelClass_HA2(name,par=par)
            model.load_est_par(f'HANK_est_par_best')
            model.load_x0s(f'HANK_2A_x0s')
        else:
            raise Exception("Wrong HH type chosen!")

    model.load_data()
    
    if ss: model.find_ss()

    if jac:
        skip_hh = skip_hh or (len(model.outputs_hh)==0)
        model.compute_jacs(skip_shocks=skip_shocks,skip_hh=skip_hh)
    
    if find_transition_path:
        model.find_transition_path_foreign_shock()

    return model

def copies(basemodel,pardicts,do_print=True,find_ss=False,skip_hh=True,compress_full=True):

    models = [basemodel]
    for i,pardict in enumerate(pardicts):

        if do_print: print(pardict)
        
        try:

            model = basemodel.copy()

            model.set_from_est_par(pardict)

            if find_ss: model.find_ss()
            model.compute_jacs(skip_shocks=True,skip_hh=skip_hh)
            model.find_transition_path_foreign_shock()

            if compress_full: model.compress_full()
            models.append(model)   
            
        except Exception as e:

            print(e)
            
    return models    

class HANKModelClass(EconModelClass, GEModelClass):

    #########
    # setup #
    #########      
     
    def set_HH_type(self,):
        """ set household type """

        self.grids_hh = ['a']  # grids
        self.pols_hh = ['a']  # policy functions
        self.inputs_hh_z = []  # transition matrix inputs
        self.inputs_hh = ['ra', 'wnT', 'wnNT', 'UniformT', 'tau', 'eps_beta', 'd', 'AI_lag']  # inputs to household problem
        self.intertemps_hh = ['vbeg_a']  # inputs to household problem
        self.outputs_hh = ['a', 'c', 'c_T', 'c_NT', 'inc_T', 'inc_NT', 'uc_T', 'uc_NT']  # output of household problem
        self.solve_hh_backwards = household_problem.solve_hh_backwards_1A

    def settings(self):
        """ fundamental settings """

        # a. namespaces
        self.namespaces = ['par', 'sol', 'sim', 'ss', 'path', 'ini']

        # b. other attributes (to save them)
        self.other_attrs = ['data']

        # c. households
        self.set_HH_type()
   
        # d. GE
        self.shocks = ['Y_s','C_s','piF_s','iF_s',
                       'G_exo','eps_beta','UniformT_exo','di','FD','ND','VAT','subP','subI']  # exogenous inputs
        
        self.unknowns = ['pi','ZT','ZNT','NFA','B','piF','piH_s','ITp','INTp','mT','mNT','wT','wNT']  # endogenous inputs

        self.targets = ['goods_mkt_T','goods_mkt_NT','NKWPCT',
                        'NKWPCNT','NFA_target',
                        'Taylor','G_budget','IT_FOC','INT_FOC','mT_res','mNT_res','Lfoc_T','Lfoc_NT']  # targets

        self.blocks = [ # list of strings to block-functions
            'blocks.pricing',
            'blocks.production_block',
            'blocks.asset_returns',
            'blocks.public_finances',
            'blocks.simple_HHs',
            'blocks.mon_pol',
            'hh', # household block
            'blocks.block_post',
            'blocks.misc_block',
            ]
        
        # Estimation
        est = self.est = SimpleNamespace()

        est.target_varlist = ['YT','YNT','I','N','labor_comp', 'Exports','r_ann','pi_ann','Q', 'piH_ann']
        est.target_T = 12

        est.subjective_weights = { # undocumented
            'C' : 1, 
            'I' : 1,  
            'N' : {'0' : 0.01}, 
            'r_ann' : {'0' : 0.1},
            'YT' : {'0' : 0.1}, 
            'YNT' : {'0' : 0.1},
            'pi_ann' : 1.,
            'piH_ann' : 1.,
            'Exports' : {'0' : 0.1},
            'CT' : {'0' : 0.01},
            'CNT' : {'0' : 0.1},
            'labor_comp' : {'0' : 0.1},
            'Q' : 1.,
        }

        self.est.pardict = {
            'NKslope__T' : (np.nan, 0.01, 1.5),
            'NKslope_NT' : (np.nan, 0.01, 1.5),
            'NKWslope__T' : (np.nan, 0.00005, 0.7),
            'NKWslope_NT' : (np.nan, 0.00005, 0.7),
            'pi_index' : (np.nan, 0.0, 0.98),
            'piW_index' : (np.nan, 0.0, 0.98),
            'phi' : (np.nan, 1.2, 2.8),
            'phi_back' : (np.nan, 0.0, 0.98),
            'phi_N' : (np.nan, 0., 12.5),
            'phi_I' : (np.nan, 0.1, 30.0),
            'kappa_r' : (np.nan, 0.0, 170),            
            'UIP_dev' : (np.nan, 0., 10.0),
            'gamma' : (np.nan, 0.1, 10.0),
            'phi_X' : (np.nan, 0., 0.99),
        }        

    def setup(self):
        """ set baseline parameters """

        par = self.par

        par.HH_type = 'HA'
        par.No_HA_list = np.array(['RA-CM', 'RA-IM', 'TA-CM', 'TA-IM'])  # RANK/TANK model labels

        par.use1A_illiquid = False # use 1A illiquid grid
        par.use_RA_jac = False  # use RA jacobian when solving

        # Preferences
        par.inv_frisch = 2.0  #* Frisch elasticity
        par.CRRA = 2.0  #* CRRA coefficient  
        
        par.beta_mean = np.nan  # discount factor, mean, range is [mean-width,mean+width] [determined in ss]
        par.beta_delta = np.nan  # discount factor, width, range is [mean-width,mean+width] [determined in ss]
        par.W2INC_target = 10.0  #* steady state target for household wealth to GDP
        par.MPC_est = np.array([0.51, 0.17, 0.11, 0.043, 0.029, 0.0267]) # From figure 4 in Fagereng et al. 2021: "MPC heterogeneity and household balance sheets"
        par.psi = np.array([np.nan,np.nan])  # disutility of work [determined in ss]

        par.alphaF = np.nan  # share of CF [determined in ss]
        par.eta = 1.5  # substitution between CF and CH
        
        par.alphaT = 0.41  #* share of CT
        par.etaT = 1.5 # substitution between CT and CNT

        # two-asset model
        par.xi = 0.0  # illiquid asset premium
        par.chi = 0.010  # parameter used to induce stationary for illiquid account. Taken From Micro jumps, Macro humps
        par.AI_target = np.nan  # Steady state level of illiquid assets [determined in ss]
        par.rai_ss_target = np.nan  # steady state level for illiquid asset return [determined in ss]

        # Earnings states and MPCs
        par.Ns = 2  # number of sectors households can be employed in
        par.Ne = 6  # number of productivity states 
        
        par.sigma_e = 0.25  #* std. of persistent shock               
        par.rho_e = 0.95  #* AR(1) parameter
        par.sT = 0.30  #* share of households that work in tradable sector

        # Mutual fund
        par.UIP_dev = np.nan  # UIP deviation [estimated]
        par.r_debt_elasticity = 0.  #* risk premium on domestic debt in RANK (closes the model) [no longer used]
        par.epsNFA = 0.5 # speed-of-adjustment
        par.tauNFA =  150 # beginning sadjustment period 
        par.deltaNFA = 30 # number of adjustment periods

        # Producer (sector-specific final good)
        par.sigma_XH = 1.0  #* substitution between capital-labor aggregate and intermediate goods
        par.alphaX = np.array([np.nan,np.nan])  # output elasticity of intermediate goods aggregate in production [determined in ss]
        par.FixedCost = np.array([np.nan,np.nan])  # fixed cost of production [determined in ss]

        par.prodbeta = np.array([np.nan,np.nan])  # aux. parameters in production function [determined in ss]
        par.prodalphaX = np.array([np.nan,np.nan])  # parameters in production function [determined in ss]

        # Intermediate goods + trade
        par.etaX = 0.5  #* substitution between different intermediate goods

        par.X_expshare = np.array([0.71, 0.46]) #* intermediate goods share
        par.K_expshare = np.array([0.10, 0.18]) #* capital share
        
        par.X_share_NT = np.array([0.60,1-0.60])  #* [NT->NT, T->NT]
        par.X_share_T = np.array([0.62,np.nan, np.nan])  # [T->T, NT->T, F->T] 

        par.HH_importshare = 0.25 #* households imports make up 25% of total imports 
        par.X2Y_target = 0.39 #* target for export/GDP ratio

        # Capital-labor
        par.use_capital = True
        par.TFP = np.nan  # aggregate productivity [determined in ss]
        par.TFP_s = np.array([np.nan,np.nan])  # sector specific productivity [determined in ss]
        par.sigma_NK = 0.5 #* substituion between employment and capital 
        par.alphaK = np.array([np.nan,np.nan])  # output elasticity of capital aggregate in production [determined in ss]

        par.prodbetaN = np.array([np.nan,np.nan])  # parameters in production function [determined in ss]
        par.prodalphaK = np.array([np.nan,np.nan])  # parameters in production function [determined in ss]

        # NKPC
        par.mu = 1.2*np.ones(2)  # price markups
        par.epsilon = np.array([np.nan,np.nan]) # substitution determining price markup [determined from mu]

        par.NKslope = np.array([np.nan,np.nan]) # slope of NKPC [estimated]
        par.theta = np.array([np.nan,np.nan])  # Rotemberg parameter [determined from NKslope+epsilon]

        par.pi_index = np.nan  # price indexation [estimated]

        # Labor firms
        par.phi_N = np.nan # adjustment cost on labor [estimated]        

        # Capital firms
        par.phi_I = np.nan # investment adjustment cost
        par.delta_K = 0.013 #* depreciation rate [Micro Jumps, Macro Humps]
        
        par.mu_r = 1.0 # outside risk premium effect
        par.kappa_r = np.nan # inner risk premium effect [estimated]

        # NKPW
        par.muw = 1.2*np.ones(2)  # markups
        par.epsilon_w = np.array([np.nan,np.nan])  # substitution determining wage markup [determined from muw]

        par.NKWslope = np.array([np.nan,np.nan]) # slope of NKPW [estimated]
        par.theta_w = np.array([np.nan,np.nan])  # Rotemberg parameter [determined from NKWslope+epsilon_w]

        par.piW_index = 0.0 # wage indexation

        # Foreign
        par.iF_s_exo = 0.005  # steady state real interest rate      

        par.alpha_s = np.nan  # scaling of foreign demand (ss.C_ast = 1.0) [determined in ss]
        par.gamma = np.nan #* Armington elasticity [estimated]
        par.phi_X = np.nan # rigidity in Armington [estimated]

        # Monetary policy
        par.MonPol = 'Taylor'  # Monetary policy rule
        par.TaylorType = 'ppi'  # Target in Taylor rule
        par.floating = True  # Floating or fixed exchange rate
        par.phi = np.nan  # Taylor rule coefficient [estimated]
        par.phi_back = np.nan  # Degree of interest rate smoothing in Taylor rule [estimated]

        # Public sector
        par.G_GDP_ratio = 0.17  #* G/GDP
        par.B_GDP_ratio = 0.95 * 4  #* B/GDP (annual)

        par.sGT_ss = 0.20  # share of public consumption going to tradeables in steady state
        par.sGT = 0.20  # share of public consumption going to tradeables when shocking G_eps

        par.sTransferT = 0.5  # share of transfers going to households in tradeables

        # Debt rule
        par.debt_rule = False  #* use public debt rule from Auclert et al 2021
        par.epsB = 1.0 #* speed-of-adjustment
        par.tauB = 50  #* adjustment period
        par.deltaB = 20 #* no adjustment periods 
        
        # Fiscal devaluation
        par.FD_shock = False  # true: VAT and subP depend on FD
        par.VAT_weight = np.nan  # weight on VAT in FD
        par.subI_weight = np.nan # subsidy on investment in FD
        par.subP_weight = np.nan # subsidy on production in FD
        
        # Grids
        self.set_HH_grids()

        # Shocks
        par.scale = np.nan  # used to scale shocks for plots

        # Misc.
        par.T = 300  # length of path
        
        par.max_iter_solve = 50_000  # maximum number of iterations when solving
        par.max_iter_simulate = 50_000  # maximum number of iterations when simulating
        par.max_iter_broyden = 100  # maximum number of iteration when solving eq. system

        par.tol_solve = 1e-11  # tolerance when solving
        par.tol_simulate = 1e-12  # tolerance when simulating
        par.tol_broyden = 1e-07  # tolerance when solving eq. system

        par.py_hh = False   
        par.py_blocks = False

        # Initial values for steady state solvers
        par.x0 = np.array([0.37189705, 0.2821916 , 1.36659349, 1.9869943 , 0.13119899])

    def set_HH_grids(self):
        """ set grids for household problem """

        par = self.par

        par.Nbeta = 3  # discount factor, number of states

        par.a_min = 0.0  # lower limit
        par.a_max = 2000  # upper limit
        par.Na = 300  # number of grid points

        par.ai_min = np.nan  # lower limit
        par.ai_max = np.nan  # upper limi
        par.Nai = 0  # number of grid points

        par.x0_het = np.array([0.9750177935216898, 0.019495793354351063]) # initial values for steady state calibration

    def allocate(self):
        """ allocate model """

        par = self.par
        par.Nfix = par.Nbeta  # number of fixed states
        par.Nz = par.Ne * par.Ns  # total number of earnings states
        
        allowed_HH_types = ['HA','HA-2A', 'RA-CM', 'RA-IM', 'TA-CM', 'TA-IM']
        assert par.HH_type in allowed_HH_types, f'HH_type = {par.HH_type} not in allowed list: {allowed_HH_types}'

        # a. containers for RANK Jacobians
        par.M_Y = np.zeros((par.T, par.T))
        par.M_R = np.zeros((par.T, par.T))
        par.M_beta = np.zeros((par.T, par.T))
           
        # b. grids
        par.a_grid = np.zeros(par.Na)
        par.beta_grid = np.zeros(par.Nbeta)
        
        par.e_grid_ss = np.zeros(par.Ne)
        par.e_trans_ss = np.zeros([par.Ne, par.Ne])
        par.e_ergodic_ss = np.zeros(par.Ne)
        
        par.z_ergodic_ss = np.zeros([par.Ns*par.Ne])
        par.z_grid_ss = np.zeros([par.Ns*par.Ne])
        
        par.s_trans_ss = np.zeros([par.Ns, par.Ns])
        par.s_ergodic_ss = np.zeros(par.Ns)
        par.s_set = np.zeros(par.Nz)
        par.sT_vec = np.array([np.nan,np.nan])

        par.tau_grid = np.zeros(par.Ns*par.Ne)

        # d. allocate GE
        self.use_FD_shock(par.FD_shock)
        self.allocate_GE()

    prepare_hh_ss = household_problem.prepare_hh_ss

    def use_FD_shock(self,FD_shock=True,update=False):
        """ Use FD shock (or not) """

        if FD_shock:
            self.par.FD_shock = True
            shocks = [x for x in self.shocks if x not in ['VAT','subP','subI','FD']]
            shocks += ['FD']
        else:
            self.par.FD_shock = False
            shocks = [x for x in self.shocks if x not in ['VAT','subP','subI','FD']]
            shocks += ['VAT','subP','subI']

        self.shocks = shocks
        if update: self.update_aggregate_settings(shocks=shocks)
        
    def use_capital(self,update=False):
        """ Use capital in the model """

        if self.par.use_capital:

            u_list_I = ['ITp', 'INTp', 'mT', 'mNT']

            K_targets = ['IT_FOC', 'INT_FOC', 'mT_res', 'mNT_res']

            if self.par.use_capital:
                u_list = u_list_I

                unknowns = [x for x in self.unknowns if x not in u_list_I]
                targets = [x for x in self.targets if x not in K_targets]
                unknowns += u_list
                targets += K_targets

            else:

                unknowns = [x for x in self.unknowns if x not in u_list_I]
                targets = [x for x in self.targets if x not in K_targets]

            self.unknowns = unknowns
            self.targets = targets
            if update: self.update_aggregate_settings(unknowns=unknowns,targets=targets)

    ################
    # steady state #
    ################
    
    def find_ss(self, do_print=False, print_timeings=False):
        """ Find steady state """
        
        steady_state.calibrate_ss(self, do_print, print_timeings)

    def ann_MPCs(self):
        """ calculate MPCs """

        par = self.par
        ss = self.ss

        # a. shock
        self._set_inputs_hh_all_ss()
        
        # 5 percent of labor income
        self.path.UniformT[:,:] = 0.0 
        self.path.UniformT[0,0] = 0.05 * (self.ss.WT * self.ss.NT + self.ss.WNT * self.ss.NNT)
        
        dI = self.path.UniformT[0,0]      

        # b. solve and simulate    
        self.solve_hh_path()
        self.simulate_hh_path()
        
        # c. calculate MPCs
        MPC = np.zeros(self.par.T)*np.nan

        C_ss = np.sum(self.ss.c*self.ss.D)    
        for t in range(self.par.T):
            MPC[t] = (np.sum(self.path.c[t]*self.path.D[t])-C_ss)/dI
            
        ann_MPC = np.zeros(round(self.par.T/4))

        for j in range(round(self.par.T/4)):
            ann_MPC[j] = np.sum(MPC[j*4:(1+j)*4])  
        
        # d. reset
        self._set_inputs_hh_all_ss()
        
        return ann_MPC, MPC
    
    ###################
    # transition path #
    ###################
    
    def scaleshock(self,varname,size=0.01,abs=False,cumeffect=False,Nq=8,log_variance=False,sign=1):

        if cumeffect:
        
            InitShock = np.sum(getattr(self.path,varname)[:Nq] / getattr(self.ss,varname) -1)
            self.par.scale = size / InitShock     
        
        elif log_variance:

            dLog = np.log(getattr(self.path,varname)) - np.log(getattr(self.ss,varname))
            InitShock_var = np.sum(dLog**2)
            self.par.scale = sign*np.sqrt(size / InitShock_var) 

            test = (np.log(getattr(self.path,varname)) - np.log(getattr(self.ss,varname)))* self.par.scale
            assert np.isclose(np.sum(test**2), size)  

        elif abs:

            InitShock = getattr(self.path,varname)[0,0] - getattr(self.ss,varname) 
            self.par.scale = size / InitShock    

        else:

            InitShock = getattr(self.path,varname)[0,0] / getattr(self.ss,varname) -1
            self.par.scale = size / InitShock        
     
    def get_foreign_shocks_from_LP(self,scale=1.0/1000):
        """ Set foreign shock from LP IRFs """

        self.par.scale = 1.0/scale

        shocks = {}
        for shockname in ['Y_s','C_s', 'piF_s', 'iF_s']:      
            
            dX = self.data.LP_IRFS_large_smoothed[shockname][:self.par.T]
            if shockname in ['Y_s','C_s']:
                shocks[f'd{shockname}'] = getattr(self.ss,shockname)*dX*scale
            else:
                shocks[f'd{shockname}'] = dX*scale
            
        return shocks

    def find_transition_path_foreign_shock(self,scale=1.0/1000,do_print=False,do_end_check=False):
        """ Find transition path with foreign shock """

        shocks = self.get_foreign_shocks_from_LP(scale)
        self.find_transition_path(shocks=shocks,do_print=do_print,do_end_check=do_end_check)

    def save_x0s(self,filename):
        """ Save x0s """

        par = self.par

        x0s = {'x0':par.x0,'x0_het':par.x0_het}
        with open(f'saved/{filename}.pkl', 'wb') as fp:
            pickle.dump(x0s,fp)

    def load_x0s(self,filename,do_print=False,skip=None):
        """ Load x0s """

        with open(f'saved/{filename}.pkl', 'rb') as fp:
            x0s = pickle.load(fp)
                                 
        for k,v in x0s.items():
            self.par.__dict__[k] = v

    ##############
    # estimation #
    ##############

    def load_data(self):
        """ Load data """

        data = self.data = SimpleNamespace()
        T = self.par.T
        
        # a. settings
        sigma = 0.9

        data_specs = {
            'Y_star':('Y_s',7),
            'P_star':('piF_s',8),
            'R_star':('iF_s',12),
            'IM_star':('C_s',6),
            'Ysum':('Y',12),
            'Csum':('C',12),
            'I':('I',12),
            'EX':('Exports',12),
            'pi':('pi_ann',12),
            'YNT':('YNT',12),
            'CNT':('CNT',12),
            'N':('N',12),
            'IM':('Imports',12),
            'RR':('r_ann',12),
            'YT':('YT',12),
            'CT':('CT',12),
            'wN':('labor_comp',12),
            'PH':('piH_ann',12),
            'Q':('Q',12) 
        }

        # b. conversion functions
        def convert_a2q(x_a):
            x_q = (1+x_a)**(1/4) - 1
            return x_q 

        def convert_a2a(x_q):
            x_a = (1+x_q)**4 - 1  
            return x_a 

        convert_q = ['piF_s', 'iF_s']

        # c. load data
        with open('LP_IRFs.pickle', 'rb') as handle:
            LP_IRFs_ = pickle.load(handle)

        with open('LP_IRFs_sum.pickle', 'rb') as handle:
            LP_IRFs_sum = pickle.load(handle)

        # d. structure and smooth data
        LP_IRFS_large = data.LP_IRFS_large = {'IRF' : {}, 'LO' : {}, 'HI' : {}, 'SE' : {}}
        LP_IRFS_SOE = data.LP_IRFS_SOE = {'IRF' : {}, 'LO' : {}, 'HI' : {}, 'SE' : {}}
        LP_IRFS_large_smoothed = data.LP_IRFS_large_smoothed = {}
        LP_IRFS_SOE_smoothed = data.LP_IRFS_SOE_smoothed = {}

        for dataname,(modelname,zero_start) in data_specs.items():

            if 'sum' in dataname:
                LP_IRFs = LP_IRFs_sum
                sign = 1
            else:
                LP_IRFs = LP_IRFs_
                sign = -1
            
            if 'star' in dataname:
                assert '_s' in modelname
                out = LP_IRFS_large
                out_smoothed = LP_IRFS_large_smoothed
            else:
                out = LP_IRFS_SOE
                out_smoothed = LP_IRFS_SOE_smoothed

            if modelname in convert_q:
                out['IRF'][modelname] = convert_a2q(-sign*LP_IRFs[dataname]['IRF']/100)
                out['LO'][modelname] = np.vstack([
                    convert_a2q(-sign*LP_IRFs[dataname]['LO']['90%']/100),
                    convert_a2q(-sign*LP_IRFs[dataname]['LO']['68%']/100)]
                    ).T
                out['HI'][modelname] = np.vstack([
                    convert_a2q(-sign*LP_IRFs[dataname]['HI']['90%']/100),
                    convert_a2q(-sign*LP_IRFs[dataname]['HI']['68%']/100)]
                    ).T
                out['SE'][modelname] = convert_a2q(-sign*LP_IRFs[dataname]['SE']/100)
            else:
                out['IRF'][modelname] = -sign*LP_IRFs[dataname]['IRF']/100
                out['LO'][modelname] = np.vstack([
                    -sign*LP_IRFs[dataname]['LO']['90%']/100,
                    -sign*LP_IRFs[dataname]['LO']['68%']/100]
                    ).T
                out['HI'][modelname] = np.vstack([
                    -sign*LP_IRFs[dataname]['HI']['90%']/100,
                    -sign*LP_IRFs[dataname]['HI']['68%']/100]
                    ).T
                out['SE'][modelname] = -sign*LP_IRFs[dataname]['SE']/100  

            shock_extended = np.hstack([out['IRF'][modelname][:zero_start],np.zeros(5)])
            smooth_irf = ndimage.gaussian_filter1d(shock_extended,sigma)

            out_smoothed[modelname] = np.hstack([smooth_irf,np.zeros(T-smooth_irf.size)])

    def prepare_estimation(self,pardict=None,target_varlist=None,target_T=None,subjective_weights=None,max_iter=200,do_print=False):
        """ Prepare estimation """

        # a. settings
        est = self.est

        est.iter = 0
        est.best_obj_val = np.inf
        est.best_x = None

        if not pardict is None: est.pardict = pardict
        if not target_varlist is None: est.target_varlist = target_varlist
        if not target_T is None: est.target_T = target_T
        if not subjective_weights is None: est.subjective_weights = subjective_weights
        if not max_iter is None: est.max_iter = max_iter

        # b. weight matrix
        self.set_weight_matrix(subjective_weights=est.subjective_weights,do_print=do_print)

    def save_est_par(self,filename):
        """ Save estimated parameters """

        par = self.par

        x = self.get_est_x(self.est.pardict)
        est_par = {parname:x[i] for i,parname in enumerate(self.est.pardict.keys())}

        with open(f'saved/{filename}.pkl', 'wb') as fp:
            pickle.dump(est_par, fp)

    def get_est_x(self,pardict):
        """ Get estimated parameters as vector """

        x = np.nan*np.ones(len(pardict))
        for i,parname in enumerate(pardict.keys()):

            if parname in ['NKslope__T', 'NKslope_NT', 'NKWslope__T', 'NKWslope_NT']:

                parname_str = parname[:-3]
                sector_str = parname[-3:]
                
                if sector_str == '__T':
                    x[i] = self.par.__dict__[parname_str][0]
                elif sector_str == '_NT':
                    x[i] = self.par.__dict__[parname_str][1]
                else:
                    raise Exception("Must be either T or NT")

            else:

                x[i] = getattr(self.par,parname)

        return x
                
    def set_from_est_par(self,est_par):
        """ Set estimated parameters """
        
        for parname,parvalue in est_par.items():

            if parname in ['NKslope__T', 'NKslope_NT', 'NKWslope__T', 'NKWslope_NT']:
                
                parname_str = parname[:-3]
                sector_str = parname[-3:]
                if sector_str == '__T':
                    self.par.__dict__[parname_str][0] = parvalue
                elif sector_str == '_NT':
                    self.par.__dict__[parname_str][1] = parvalue
                else:
                    raise Exception("Must be either T or NT")
                
            else:

                setattr(self.par,parname,parvalue)

    def load_est_par(self,filename,do_print=False,skip=None):
        """ Load estimated parameters from saved file """

        # a. load
        with open(f'saved/{filename}.pkl', 'rb') as fp:
            est_par = pickle.load(fp)

        if skip is not None:
            for parname in skip:
                del est_par[parname]

        # b. set
        if do_print:
            
            # hack to get old values
            pardict = {parname:None for parname in est_par.keys()}
            x = self.get_est_x(pardict)

            # print
            for (parname,parvalue),oldvalue in zip(est_par.items(),x):
                print(f'{parname} = {parvalue:.4f} [now: {oldvalue:.4f}]')
                                 
        self.set_from_est_par(est_par)

    def set_weight_matrix(self,subjective_weights=None,do_print=False):
        """ Set weight matrix for estimation """

        est = self.est
        data = self.data

        # a. diag
        Sigma_diag = None
        for j,var in enumerate(est.target_varlist):
            if Sigma_diag is None:
                Sigma_diag = data.LP_IRFS_SOE['SE'][var][:est.target_T]**2
            else:
                Sigma_diag = np.hstack((Sigma_diag,data.LP_IRFS_SOE['SE'][var][:est.target_T]**2))

        Sigma = np.eye(len(est.target_varlist)*est.target_T)
        np.fill_diagonal(Sigma, Sigma_diag)
        Sigma_inv = np.linalg.inv(Sigma)

        # b. subjective weights (should perhaps be specified further out)
        if subjective_weights is None: subjective_weights = self.subjective_weights

        subj_weightmat = np.eye((len(est.target_varlist)*est.target_T))
        for j,key in enumerate(est.target_varlist):
            if key in subjective_weights:
                if isinstance(subjective_weights[key], dict):
                    for key2 in subjective_weights[key]:
                        subj_weightmat[j*est.target_T+int(key2),j*est.target_T+int(key2)] = subjective_weights[key][key2]
                else:
                    for t in range(est.target_T):
                        subj_weightmat[j*est.target_T+t,j*est.target_T+t] = subjective_weights[key]

        # c. total weight matrix
        est.Sigma_inv = Sigma_inv @ subj_weightmat
        est.var_data = Sigma

        # d. show
        if do_print:

            Sigma_inv_diag = np.diag(est.Sigma_inv)
            temp = {}
            temp['t'] = np.arange(est.target_T)
            for j,key in enumerate(est.target_varlist):
                temp[key] = Sigma_inv_diag[j*est.target_T:(1+j)*est.target_T] / 10000 # scaling of objective function 

            df = pd.DataFrame(temp)
            display(df.head(est.target_T))

    def est_obj(self,x,do_print=False,in_estimation=False):
        """ Objective function for estimation """

        t0 = time.time()

        par = self.par
        ss = self.ss
        path = self.path
        est = self.est
        data = self.data

        if in_estimation: 

            if os.path.isfile(f'est.stop'): 
                os.remove(f'est.stop')
                raise StopIteration

            if est.iter >= est.max_iter: raise StopIteration        
            est.iter += 1

        # a. set parameters
        if not x is None:

            if do_print:

                np.set_printoptions(suppress=True)
                print(f'x: {np.round(x,5)}',end='')
                np.set_printoptions(suppress=False)

            est_par = {parname:x[i] for i,parname in enumerate(est.pardict.keys())}
            self.set_from_est_par(est_par)

        # b. IRF deviations
        try:
            
            # i. evaluate
            self.compute_jacs(do_print=False,skip_shocks=False,skip_hh=True)
            self.find_transition_path_foreign_shock(do_print=False)

            # ii. deviations
            IRF_devs = np.zeros(len(est.target_varlist)*est.target_T) * np.nan 
            model_IRFs = np.zeros(len(est.target_varlist)*est.target_T) * np.nan 
            for j,var in enumerate(est.target_varlist):

                if var in ['r_ann','pi_ann', 'piH_ann']:
                    model_IRFs[j*est.target_T:(1+j)*est.target_T]  = (getattr(path, var)[:,0]-getattr(ss, var))[:est.target_T]*par.scale
                else:
                    model_IRFs[j*est.target_T:(1+j)*est.target_T]  = (getattr(path, var)[:,0]-getattr(ss, var))[:est.target_T]*par.scale/getattr(ss, var)

                IRF_devs[j*est.target_T:(1+j)*est.target_T] = model_IRFs[j*est.target_T:(1+j)*est.target_T]  - data.LP_IRFS_SOE['IRF'][var][:est.target_T]
                
            obj_val = (IRF_devs.T @ est.Sigma_inv @ IRF_devs)
            
            if do_print: print(f' -> {obj_val:.4f}',end='')

        except Exception as e:
                
            print(e)
            obj_val = np.inf

        # c. return
        if in_estimation:

            if obj_val < est.best_obj_val:
                
                est.best_obj_val = obj_val
                est.best_x = x
                self.save_est_par(filename=f'{self.name}_est_par_best_yet')
                if do_print: print(' [best yet]',end='')

        if do_print: print(f' [took {time.time()-t0:.1f} secs]')

        if in_estimation:
            return obj_val
        else:
            return obj_val, model_IRFs
    
    def est_se(self,do_print=False):
        """ Compute standard errors """

        est = self.est

        # a. baseline
        h = 1e-4
        x = self.get_est_x(self.est.pardict)
        _,model_IRFs_base = self.est_obj(x) # baseline

        # b. G
        G = np.zeros((model_IRFs_base.size,len(est.pardict)))
        for i, parname in enumerate(est.pardict):

            if do_print: print(f' {parname}')
            
            xh = deepcopy(x)
            xh[i] += h

            _,model_IRFs_alt = self.est_obj(xh) # derivative

            G[:,i] = (model_IRFs_alt - model_IRFs_base) / h 

        # c. reset
        _,model_IRFs_base = self.est_obj(x) # reverse to baseline

        # d. formula
        W = est.Sigma_inv 
        GW  = np.transpose(G) @ W
        GWG = GW @ G
        GWG_inv = np.linalg.inv(GWG)

        V = est.var_data
        Avar = GWG_inv @ (GW @ V @ np.transpose(GW) ) @ GWG_inv

        # e. save
        est.se = {}
        for i,parname in enumerate(est.pardict.keys()):
            est.se[parname] = np.sqrt(Avar[i,i])

    def estimate(self,do_se=False,do_print=True):
        """ Estimate model with impulse responses matching """

        if os.path.isfile(f'est.stop'): os.remove(f'est.stop')
        if os.path.isfile(f'est.running'): os.remove(f'est.running')

        est = self.est

        est.iter = 0
        est.best_obj_val = np.inf
        est.best_x = None

        # a. intial values and bounds
        x0 = np.nan*np.ones(len(self.est.pardict))
        bounds = [None]*len(self.est.pardict)

        x0_ = self.get_est_x(self.est.pardict)
        for i,(parname,(value,lb,ub)) in enumerate(est.pardict.items()):
            x0[i] = value if not np.isnan(value) else x0_[i]
            bounds[i] = (lb,ub)
            
        # b. run
        obj = lambda x: self.est_obj(x,do_print,True)

        try:
            open('est.running', mode='w').close()
            res = optimize.minimize(obj,x0,method='Nelder-Mead',tol=1e-04,bounds=bounds) 
        except StopIteration:
            print(f'estimation stopped after {est.iter} iterations')
        
        if os.path.isfile(f'est.running'): os.remove(f'est.running')

        if not est.best_x is None:

            # c. finalize
            print('\nfinal evaluation:')
            self.est_obj(est.best_x,do_print)

        # c. standard errors
        if do_se: 

            if do_print: 
                print('\ncalculating standard errors:')

            self.est_se(do_print=do_print)

        # d. finalize
        if do_print:

            print('')
            if not est.best_x is None: print(f'best obj val: {est.best_obj_val:.4f}')

            x = self.get_est_x(self.est.pardict)
            for i,parname in enumerate(est.pardict.keys()):

                print(f' {parname} = {x[i]:.4f}',end='')
                if do_se:
                    print(f' [s.e. {est.se[parname]:.4f}]')
                else:
                    print('')

class HANKModelClass_RA(HANKModelClass):

    def set_HH_type(self):

        self.pols_hh = []  # policy functions 
        self.grids_hh = []  # grids
        self.inputs_hh_z = []  # transition matrix inputs
        self.inputs_hh = []  # inputs to household problem
        self.intertemps_hh = []  # inputs to household problem
        self.outputs_hh = []  # output of household problem
        self.solve_hh_backwards = None          

    def set_HH_grids(self):

        par = self.par

        par.Nbeta = 1
        
        par.a_min = np.nan
        par.a_max = np.nan
        par.Na = 0
        
        par.ai_min = np.nan
        par.ai_max = np.nan
        par.Nai = 0
        
        par.x0_het = np.nan*np.ones(2)
    
class HANKModelClass_HA2(HANKModelClass):

    def set_HH_type(self):

        self.grids_hh = ['a', 'ai']  # grids
        self.pols_hh = ['a', 'ai']  # policy functions
        self.inputs_hh_z = []  # transition matrix inputs
        self.inputs_hh = ['ra', 'rai', 'wnT', 'wnNT', 'UniformT', 'tau', 'eps_beta']  # inputs to household problem
        self.intertemps_hh = ['vbeg_a']  # inputs to household problem
        self.outputs_hh = ['a', 'ai', 'd_payout', 'c', 'c_T', 'c_NT', 'inc_T', 'inc_NT', 'uc_T', 'uc_NT']  # output of household problem
        self.solve_hh_backwards = household_problem.solve_hh_backwards_2A

    def set_HH_grids(self):

        par = self.par

        par.Nbeta = 1
        
        par.a_min = 0.0
        par.a_max = 30
        par.Na = 75
        
        par.ai_min = 3.0
        par.ai_max = 8.0
        par.Nai = 20

        par.x0_het = np.array([0.9802443299010444, 5.016387077146562]) # initial values for steady state calibration      