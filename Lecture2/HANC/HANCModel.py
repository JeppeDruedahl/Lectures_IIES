import numpy as np

from EconModel import EconModelClass
from GEModelTools import GEModelClass

import steady_state
import household_problem
import blocks

class HANCModelClass(EconModelClass,GEModelClass):    

    def settings(self):
        """ fundamental settings """

        # a. namespaces (typically not changed)
        self.namespaces = ['par','ini','sim','ss','path']
        
        # b. household
        self.grids_hh = ['a'] # grids
        self.pols_hh = ['a'] # policy functions
        self.inputs_hh = ['r','w'] # direct inputs
        self.inputs_hh_z = [] # transition matrix inputs
        self.outputs_hh = ['a','c','l'] # outputs
        self.intertemps_hh = ['vbeg_a'] # intertemporal variables

        # c. GE
        self.shocks = ['Gamma'] # exogenous shocks
        self.unknowns = ['K','L'] # endogenous unknowns
        self.targets = ['clearing_A','clearing_L'] # targets = 0
        self.blocks = [ # list of strings to block-functions
            'blocks.production_firm',
            'blocks.mutual_fund',
            'hh', # household block
            'blocks.market_clearing']

        # d. functions
        self.solve_hh_backwards = household_problem.solve_hh_backwards

    def setup(self):
        """ set baseline parameters """

        par = self.par

        par.Nfix = 3 # number of fixed discrete states (here discount factor)
        par.Nz = 7 # number of stochastic discrete states (here productivity)

        # a. preferences
        par.sigma = 2.0 # CRRA coefficient
        par.beta_mean = 0.975 # discount factor, mean, range is [mean-width,mean+width]
        par.beta_delta = 0.010 # discount factor, width, range is [mean-width,mean+width]

        # b. income parameters
        par.rho_z = 0.95 # AR(1) parameter
        par.sigma_psi = 0.30*(1.0-par.rho_z**2.0)**0.5 # std. of persistent shock

        # c. production and investment
        par.alpha = 0.36 # cobb-douglas
        par.delta = np.nan # depreciation [determined in ss]

        # d. calibration
        par.r_ss_target = 0.01 # target for real interest rate
        par.w_ss_target = 1.0 # target for real wage

        # f. grids         
        par.a_max = 500.0 # maximum point in grid for a
        par.Na = 300 # number of grid points

        # g. shocks
        par.jump_Gamma = -0.10 # initial jump
        par.rho_Gamma = 0.90 # AR(1) coefficient
        par.std_Gamma = 0.01 # std. of innovation

        # h. misc.
        par.T = 500 # length of transition path        
        par.simT = 2_000 # length of simulation 
        
        par.max_iter_solve = 50_000 # maximum number of iterations when solving household problem
        par.max_iter_simulate = 50_000 # maximum number of iterations when simulating household problem
        par.max_iter_broyden = 100 # maximum number of iteration when solving eq. system
        
        par.tol_solve = 1e-12 # tolerance when solving household problem
        par.tol_simulate = 1e-12 # tolerance when simulating household problem
        par.tol_broyden = 1e-10 # tolerance when solving eq. system
        
        par.py_hh = False # call solve_hh_backwards in Python-model
        par.py_block = True # call blocks in Python-model
        par.full_z_trans = False # let z_trans vary over endogenous states

    def allocate(self):
        """ allocate model """

        par = self.par

        # a. grids
        par.Nbeta = par.Nfix
        par.beta_grid = np.zeros(par.Nbeta)
        
        # b. solution
        self.allocate_GE() # should always be called here

    prepare_hh_ss = steady_state.prepare_hh_ss
    find_ss = steady_state.find_ss