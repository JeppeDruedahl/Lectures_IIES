import numpy as np

import matplotlib.pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

from EconModel import EconModelClass
from GEModelTools import GEModelClass

import household_problem
import steady_state

class HANKSAMModelClass(EconModelClass,GEModelClass):    

    #########
    # setup #
    #########

    def settings(self):
        """ fundamental settings """

        # a. namespaces
        self.namespaces = ['par','ini','ss','path','sim']
        
        # b. household
        self.grids_hh = ['a']
        self.pols_hh = ['a']
        self.inputs_hh = ['w','RealR_ex_post','tau']
        self.inputs_hh_z = ['delta','lambda_u']
        self.outputs_hh = ['a','c','s']
        self.intertemps_hh = ['vbeg_a']

        # c. GE
        self.shocks = ['shock_TFP','delta','w']
        self.unknowns = ['px','Vj','vt','ut','S','Pi']

        self.targets = ['errors_Vj','errors_ut','errors_entry','errors_Pi','errors_assets','errors_search']

        # d. functions
        self.solve_hh_backwards = household_problem.solve_hh_backwards
        self.blocks = [
            'blocks.production',
            'blocks.labor_market',
            'blocks.entry',
            'blocks.price_setters',
            'blocks.central_bank',
            'blocks.government',
            'blocks.mutual_fund',
            'hh',
            'blocks.market_clearing']

    def setup(self):
        """ set baseline parameters """

        par = self.par
        par.Nfix = 3

        # a. targets
        par.RealR_ss = 1.02**(1/12)
        par.mutual_fund_share = 0.0
        par.div_tax_share_ss = 0.90
        par.div_tax = np.nan

        # b. preferences
        par.beta_HtM = 0.92**(1/12)
        par.beta_mid = 0.94**(1/12)
        par.beta_PIH = 0.975**(1/12)

        par.HtM_share = 0.30
        par.PIH_share = 0.30

        par.beta = 0.96**(1/12) # discount factor
        par.sigma = 2.0 # CRRA coefficient

        # b. matching and bargaining
        par.A = np.nan # matching efficiency, determined endogenously
        par.theta_ss = 0.60 # tightness in ss
        par.alpha = 0.60 # matching elasticity
        par.lambda_u_ss = 0.30
        
        # c. intermediary goods firms
        par.w_share_ss = 0.97 # wage in steady state
        par.kappa = np.nan # flow vacancy cost, determined endogenously
        par.delta_ss = 0.025

        # d. final goods firms
        par.epsilon_p = 20.0 # price elasticity      
        par.phi = 600.0 # Rotemberg cost

        # e. monetary policy
        par.rho_R = 0.0 # inertia
        par.delta_pi = 1.5 # inflation aggressiveness

        # f. government
        par.qB_share_ss = 1.00 # government bonds (share of wage)

        par.UI_ratio = 0.25 # UI ratio
    
        par.omega = 0.90 # responsiveness of tax to debt
        par.delta_q = 1-1/60 # maturity of government bonds

        # g. shocks
        par.rho_shock_TFP = 0.965 # persitence
        par.jump_shock_TFP = -0.007 # jump

        # h. household problem
        par.Nz = 2 # number of income states
        par.Na = 500 # number of asset grid points
        par.a_max = 200 # max level of assets

        # d. misc
        par.T = 300 # length of path        
        
        par.max_iter_solve = 50_000 # maximum number of iterations when solving
        par.max_iter_simulate = 50_000 # maximum number of iterations when simulating
        par.max_iter_broyden = 50 # maximum number of iteration when solving eq. system
        
        par.tol_solve = 1e-12 # tolerance when solving
        par.tol_simulate = 1e-12 # tolerance when simulating
        par.tol_broyden = 1e-8 # tolerance when solving eq. system
        par.tol_R = 1e-12 # tolerance when finding RealR for ss
        par.tol_calib = 1e-5 # tolerance when calibrating (C_drop and var_u)

        par.py_hh = False
        par.py_blocks = False
        par.full_z_trans = True

    def allocate(self):
        """ allocate model """
        
        par = self.par        

        par.beta_grid = np.zeros(par.Nfix)
        par.beta_shares = np.zeros(par.Nfix)

        self.allocate_GE()

    prepare_hh_ss = steady_state.prepare_hh_ss
    find_ss = steady_state.find_ss