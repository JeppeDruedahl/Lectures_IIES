import os
import pickle
from copy import deepcopy
import time
import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

from EconModel import EconModelClass
from GEModelTools import GEModelClass
from consav.misc import elapsed

import household_problem
import steady_state
import blocks
import duration_dependent_UE

class StopIteration(Exception): pass

class FullHANKSAMModelClass(EconModelClass,GEModelClass):    

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
        self.inputs_hh = ['shock_beta','w','RealR_ex_post','tau','u_bar','phi_obar','hh_div','hh_transfer']
        self.inputs_hh_z = ['delta','lambda_u_s']
        self.outputs_hh = ['a','c','v_','s','u_UI','u_ALL']
        self.intertemps_hh = ['vbeg','vbeg_a']

        # c. GE
        self.shocks = ['shock_TFP','shock_beta','u_bar','phi_obar','retention_subsidy','hiring_subsidy','hh_transfer','G','u_target']
        self.unknowns = ['px','Vj','Vv','Pi','u','vt','U_UI_hh_guess','S']
        self.targets = ['errors_Vj','errors_Vv','errors_Pi','errors_assets','errors_u','errors_vt','errors_U_UI','errors_search']
        self.blocks = [
            'blocks.wage_and_profits',
            'blocks.Vj',
            'blocks.entry',
            'blocks.labor_market',
            'blocks.Vv',
            'blocks.wage_rule',
            'blocks.phillips',
            'blocks.central_bank',
            'blocks.dividends',
            'blocks.arbitrage',
            'blocks.government',
            'blocks.hh_vars',
            'hh',
            'blocks.market_clearing']

        # d. functions
        self.solve_hh_backwards = household_problem.solve_hh_backwards

        # e. misc
        self.other_attrs = ['data','moms','varlabels']
        self.varlabels = {
            'v':'vacancy stock',
            'u':'unemployment rate',
            'delta':'separation rate',
            'theta':'tightness',
            'lambda_u':'job-finding-rate',
            'lambda_v':'job-filling rate',
            'shock_TFP':'TFP',
            'px':'intermediary goods price',
            'R':'nominal Interest rate',
            'Pi':'inflation',
            'RealR':'real interest rate',
            'RealR_ex_post':'real interest rate (ex post)',
            'q':'price of government debt',
            'qB':'value of government debt',
            'tau':'tax rate'
        }

    def setup(self):
        """ set baseline parameters """

        par = self.par

        # a. settings
        par.Nfix = 3 # number of betas
        
        # b. macros
        par.RA = False # representative agent
        par.wage_setting = 'fixed' # wage setting
        par.free_entry = False # free entry
        par.exo_sep = False # exogenous separation
        par.do_vbeg = False # calculate beginning of period value function

        # c. targets
        par.RealR_ss = 1.02**(1/12) # real interest

        par.UEs_taget = np.array([0.0,-0.56,-0.68,-0.7,-0.72,-0.74,-0.76,-0.775,-0.8,-0.82,-0.85,-0.87,-0.90]) # log-difference of duration dependent UE rates
        par.UI_of_u_share_target = 39 # UI share of unemployed

        par.C_drop_ss_target = -20 # C of unemployed vs employed
        par.C_drop_ex_ss_target = -43 # C drop of unemployed at expiration rel. to income drop
        par.MPC_qtr_target = 40 # quarterly MPC

        par.timeshift_target = 9 # top of EU vs. bottom of UE (must be integer)
        par.var_u_target = 0.88 # variance of unemployment to 1 percent TFP shock
        par.EU_share = 45.0 # EU share of unemployment variance in decomposition

        par.eps_calib_ss = 0.01 # stop calibration in ss if below
        par.eps_calib_path = 1.0 #0.10 # stop calibration of path if below

        # d. model parameters

        # preferences
        par.sigma = 2.0 # CRRA coefficient

        par.beta_min = 0.95**(1/12) # low beta type
        par.beta_max = 0.975**(1/12) # high beta type
        
        par.HtM_share = 0.30 # share of HtM
        par.beta_HtM_cutoff = 0.90 # always HtM if below

        par.beta_PIH = (1/par.RealR_ss**12-0.005)**(1/12) # beta of PIH
        par.PIH_share = 0.0

        par.beta_firm = 0.98**(1/12) # discount factor of firms     

        par.vartheta = 0.01 # disutility of work
        par.zeta = 1.0 # disutility of search effort, scale
        par.nu = 20.0 # disutility of search effort, curvature

        # serach bounds
        par.beta_grid_min = 0.95
        par.beta_grid_max = (0.995*1/par.RealR_ss**12)**(1/12)
        
        # searching
        par.exo_search = True

        # firms
        par.adj_virtual_share = 1.0 # share of adjustment costs considered virtual
        par.div_hh = 1.0 # share of dividends distributed to households
        par.div_PIH = 0.0 # share of dividends distributed to PIH
        par.div_return  = False # dividends distributed relative to asssets
        
        # matching and bargaining
        par.A = np.nan # matching efficiency, determined endogenously
        par.theta_ss = 0.60 # tightness in ss
        par.alpha = 0.60 # matching elasticity
        par.rho_w = 0.00 # wage rule, elasticity
        par.eta_u = 0.00 # wage rule, elasticity
        par.eta_e = 0.00 # wage rule, elasticity
        par.eta_TFP = 0.00 # wage rule, elasticity

        # intermediary goods firms
        par.w_ss = 0.50 # wage in steady state
        par.kappa = np.nan # flow vacancy cost, determined endogenously
        par.kappa_0 = 0.1 # fixed vacancy cost
        par.psi = 1.0 # separation elasticity
        par.xi = 0.010 # entry elasticity
        par.p_fac = 1.20 # factor for maximum increase in separation rate
        par.p = np.nan # maximum separation rate, determined endogenously
        par.Upsilon = np.nan # Vj at maximum separation rate, determined endogenously

        # final goods firms
        par.epsilon_p = 6.0 # price elasticity      
        par.phi = 355.0 # Rotemberg cost (Calvo of 9 month)

        # monetary policy
        par.rho_R = 0.0 # inertia
        par.delta_pi = 1.5 # inflation aggressiveness

        # government
        par.Nu = 13 # number of u states
        par.NUI = 2 # number of UI states

        par.u_bar_ss = 6.0 # UI duration
        par.phi_obar_ss = 0.76 # high UI ratio (rel. to w) *before* exhausation (in steady state)
        par.phi_ubar = 0.55 # low UI ratio (rel. to w) *after* exhausation
        par.UI_prob = np.nan # probability of getting UI, determined endogenously

        par.tau_ss = 0.30 # tax rate
        par.omega = 0.10 # responsiveness of tax to debt
        par.delta_q = 1-1/60 # maturity of government bonds

        # e. shocks
        par.rho_shock_TFP = 0.907**(1/3) # persitence
        par.jump_shock_TFP = -0.01 # jump

        par.rho_shock_beta = 0.965 # persistence
        par.jump_shock_beta = (1.0+0.0108)**(1/12)-1.0 # jump

        par.rho_G = 0.965 # persistence
        par.jump_G = 0.001 # jump

        # f. household problem
        par.Na = 100 # number of asset grid points
        par.a_max = 100 # max level of assets

        # g. misc
        par.T = 600 # length of path        
        
        par.max_iter_solve = 500_000 # maximum number of iterations when solving
        par.max_iter_simulate = 500_000 # maximum number of iterations when simulating
        par.max_iter_broyden = 50 # maximum number of iteration when solving eq. system
        
        par.tol_solve = 1e-12 # tolerance when solving
        par.tol_simulate = 1e-12 # tolerance when simulating
        par.tol_broyden = 1e-8 # tolerance when solving eq. system
        
        par.py_hh = False # call solve_hh_backwards in Python-model
        par.py_block = False # call blocks in Python-model
        par.full_z_trans = True # let z_trans vary over endogenous states

    def create_beta_grid(self):
        """ create beta grid """

        par = self.par

        # a. grid
        par.beta_grid = np.zeros(par.Nfix)
        par.beta_grid[0] = 0.0

        if par.Nfix == 3:
            par.beta_grid[1:] = np.array([par.beta_max,par.beta_PIH]) 
        elif par.Nfix > 3:
            beta = par.beta_grid_min
            for i_beta in range(1,par.Nfix):
                    par.beta_grid[i_beta] = beta
                    beta = 0.5*beta + 0.5*par.beta_grid_max
        else:
            raise ValueError('Nfix must be at least 3')

        # b. shares
        par.beta_shares = np.zeros(par.Nfix)
        par.beta_shares[0] = par.HtM_share
        if par.Nfix == 3:
            par.beta_shares[1] = 1-par.HtM_share-par.PIH_share
            par.beta_shares[2] = par.PIH_share
        else:
            par.beta_shares[1:] = (1-par.HtM_share)/(par.Nfix-1)
 
    def create_grids(self):
        """ create grids """

        par = self.par

        # a. z
        par.Nz = par.Nu*par.NUI+1

        par.i_u_hh = np.zeros((par.Nfix,par.Nz))
        par.i_UI_hh = np.zeros((par.Nfix,par.Nz))

        for i_fix in range(par.Nfix):
            par.i_u_hh[i_fix,0] = 0
            par.i_UI_hh[i_fix,0] = 0
            i_z = 1
            for i_UI in range(par.NUI):
                for i_u in range(par.Nu):
                    par.i_u_hh[i_fix,i_z] = 1+i_u
                    par.i_UI_hh[i_fix,i_z] = i_UI
                    i_z += 1

        par.s_eff = np.zeros(par.Nu+1)

        # b. beta
        self.create_beta_grid()

    def allocate(self):
        """ allocate model """
        
        # a. grid
        self.create_grids()

        # c. GE
        self.allocate_GE() # should always be called here

        # d. moms
        self.moms = {}

    prepare_hh_ss = steady_state.prepare_hh_ss
    find_ss = steady_state.find_ss

    ##############
    # more setup #
    ##############

    def save(self,name=None):
        """ save model """

        if name is None: name = self.name

        with open(f'saved/par_{name}.p', 'wb') as f:
            pickle.dump(self.par.__dict__, f)

    def load_data(self):
        """ load data """

        with open('data/data.p', 'rb') as f:
            self.data = pickle.load(f)

    def load(self,name=None):
        """ load parameters from file """

        if name is None: name = self.name

        # a. load
        with open(f'saved/par_{name}.p', 'rb') as f:
            par_dict_loaded = pickle.load(f)

        # b. update parameters
        for key in par_dict_loaded.keys():
            self.par.__dict__[key] = par_dict_loaded[key]

        # c. set macros
        self.set_macros_auto()

        # d. load data
        self.load_data()
    
    def copy(self,name=None,**kwargs):
        """ copy model - with auto setting of macros """

        other = super().copy(name=name,**kwargs)

        return other

    def set_macros(self,free_entry=None,wage_setting=None):
        """ set macros """

        par = self.par

        # baseline
        unknowns = ['px','Vj','Vv','Pi','u','vt','U_UI_hh_guess','S']
        targets = ['errors_Vj','errors_Vv','errors_Pi','errors_assets','errors_u','errors_vt','errors_U_UI','errors_search']

        if not free_entry is None: par.free_entry = free_entry
        if not wage_setting is None: par.wage_setting = wage_setting 

        # a. free entry
        if par.free_entry:
            unknowns = ['px','entry','Vj','Pi','u','vt','U_UI_hh_guess','S']

        # b. wage setting
        if par.wage_setting == 'fixed':
            pass
        elif par.wage_setting == 'rule':
            unknowns += ['w']
            targets += ['errors_WageRule']
        else:
            raise NotImplementedError

        self.update_aggregate_settings(unknowns=unknowns,targets=targets)

    def set_macros_auto(self):
        """ set macros automatically """

        par = self.par
        if np.isclose(self.par.psi,0.0): self.par.exo_sep = True

        do_free_entry = np.isinf(self.par.xi)
        do_wage_rule = not np.isclose(self.par.eta_u,0.0) or not np.isclose(self.par.eta_e,0.0) or not np.isclose(self.par.eta_TFP,0.0)
            
        if do_free_entry and do_wage_rule:
            self.set_macros(free_entry=True,wage_setting='rule')
        elif do_free_entry:
            self.set_macros(free_entry=True,wage_setting='fixed')
        elif do_wage_rule:
            self.set_macros(free_entry=False,wage_setting='rule')

    def get_IRF(self,varname):
        """ get IRF from model or data """

        # a. steady state
        ssvalue = getattr(self.ss,varname)

        # b. values
        values = getattr(self.path,varname)[:,0]

        # c. transformation
        if varname in ['R','RealR','Pi']:

            IRF = 100*(values**12-ssvalue**12)
            ylabel = '%-points (ann.)'

        elif varname in ['u','delta','lambda_u','lambda_v','shock_beta']:

            IRF = 100*(values-ssvalue)
            ylabel = '%-points'

        elif varname in ['u_bar']:

            IRF = values-ssvalue
            ylabel = 'months'

        else:

            IRF = 100*(values - ssvalue)/np.abs(ssvalue)
            ylabel = '%'

        # d. return
        return IRF,ylabel
    
    def get_flex(self,do_print=False):
        """ create copy with flexible prices """

        model_flex = self.copy()
        model_flex.par.phi = 0.001
        model_flex.compute_jacs(skip_hh=True,skip_shocks=True)
        model_flex.find_transition_path(shocks=['shock_TFP'],do_end_check=False)
        model_flex.calc_moms_path(do_print=do_print)

        return model_flex
    
    def get_RA(self,do_print=False):
        """ create copy with RA """

        model_RA = self.copy()
        model_RA.par.RA = True
        model_RA.compute_jacs(skip_hh=True,skip_shocks=True)
        model_RA.find_transition_path(shocks=['shock_TFP'],do_end_check=False)
        model_RA.calc_moms_path(do_print=do_print)            

        return model_RA
    
    def full_run(self,do_ss=True,calib_beta=False,recalib_beta=False,calib_path=False,calib_path_ini=False,
                 skip_hh=False,do_print_full=False,do_print=False,do_save=False):
        """ full run of model """

        par = self.par

        if par.exo_search:
            self.calibrate_dur_dep_UE()

        if do_ss: 
            
            self.find_ss(calib_beta=calib_beta,recalib_beta=recalib_beta,do_print=do_print_full)
            if do_print_full: print('')
        
        if calib_path: 
        
            if calib_path_ini:
                x0 = self.calibrate_path_ini(do_print=do_print_full)
            else:
                x0 = np.array([par.psi,par.xi,par.w_ss])

            self.calibrate_path(x0,do_print=do_print_full)
            if do_print_full: print('')

        self.compute_jacs(skip_hh=skip_hh,skip_shocks=True)

        try:
            self.find_transition_path(shocks=['shock_TFP'],do_end_check=False)
            self.calc_moms_ss()
            self.calc_moms_path(do_print=do_print)
        except Exception as e:
            pass

        if do_print_full: print('')
        if do_save: self.save()

    ###############
    # calibration #
    ###############

    calibrate_dur_dep_UE = duration_dependent_UE.calibrate 

    def calc_Cs(self,i_fix=None):
        """ calculate consumption """

        par = self.par
        ss = self.ss

        if i_fix is None:
            i = 0
            j = par.Nfix
        else:
            i = i_fix
            j = i + 1
        
        C_u_dur = np.nan*np.ones(par.Nu)

        if np.sum(ss.D[i:j]) > 0.0:

            C_e = np.sum(ss.D[i:j,0,:]*ss.c[i:j,0,:])/np.sum(ss.D[i:j,0,:])
            C_u = np.sum(ss.D[i:j,1:,:]*ss.c[i:j,1:,:])/np.sum(ss.D[i:j,1:,:])

            for i_z in range(par.Nz):

                i_u = par.i_u_hh[i:j,i_z]
                i_UI = par.i_UI_hh[i:j,i_z]
                if all(i_UI) == 0: continue
                assert all(i_UI) == 1
                C_u_dur[int(i_u[0]-1)] = np.sum(ss.D[i:j,i_z,:]*ss.c[i:j,i_z,:])/np.sum(ss.D[i:j,i_z,:])
        
        else:

            C_e = np.nan
            C_u = np.nan

        return C_e,C_u,C_u_dur
        
    def calc_moms_ss(self,vec=False,do_print=False):
        """ calculate steady state moments """

        t0 = time.time()

        par = self.par
        ss = self.ss
        moms = self.moms

        # a. MPC
        model_MPC = self.copy()
        model_MPC.par.T = 3 # quarterly MPC
        
        model_MPC.allocate_GE()
        model_MPC.par.z_grid = self.par.z_grid # overwritten
        model_MPC.par.a_grid = self.par.a_grid # overwritten
        model_MPC.ss = deepcopy(self.ss) # overwritten
        
        hh_transfer_ss = ss.hh_transfer
        dy = 0.01*(1-ss.tau)*ss.w # 1 percent of after-tax labor income

        custom_paths = {'hh_transfer':hh_transfer_ss*np.ones(3)}
        custom_paths['hh_transfer'][0] += dy
    
        path = model_MPC.decompose_hh_path(do_print=False,use_inputs=['hh_transfer'],custom_paths=custom_paths);
        moms['MPC_qtr'] = 100*np.sum(path.C_hh[:3,0]-ss.C_hh)/dy

        if vec:

            moms['MPC_qtr_vec'] = np.zeros(par.Nfix)
            for i_fix in range(par.Nfix):
                
                if np.sum(ss.Dbeg[i_fix]) > 0.0:

                    Dbeg = ss.Dbeg.copy()
                    Dbeg[:] = 0.0
                    Dbeg[i_fix] = ss.Dbeg[i_fix]
                    Dbeg /= np.sum(Dbeg)

                    path_base = model_MPC.decompose_hh_path(do_print=False,Dbeg=Dbeg,use_inputs=[],custom_paths=custom_paths);
                    path = model_MPC.decompose_hh_path(do_print=False,Dbeg=Dbeg,use_inputs=['hh_transfer'],custom_paths=custom_paths);

                    moms['MPC_qtr_vec'][i_fix] = 100*np.sum(path.C_hh[:4,0]-path_base.C_hh[:4,0])/dy

                else:

                    moms['MPC_qtr_vec'][i_fix] = np.nan

        # b. consumption moments 
        C_e,C_u,C_u_dur = self.calc_Cs()

        moms['C_drop_ss'] = (C_u/C_e-1)*100

        C_drop_ex = (C_u_dur[6]-C_u_dur[5])/((1-ss.tau)*(ss.phi_obar-par.phi_ubar)*ss.w)*100
        moms['C_drop_ex'] = C_drop_ex

        if vec:

            moms['C_drop_ss_vec'] = np.zeros(par.Nfix)
            moms['C_drop_ex_vec'] = np.zeros(par.Nfix)

            for i_fix in range(par.Nfix):

                C_e,C_u,C_u_dur = self.calc_Cs(i_fix=i_fix)
                moms['C_drop_ss_vec'][i_fix] = (C_u/C_e-1)*100
                moms['C_drop_ex_vec'][i_fix] = (C_u_dur[6]-C_u_dur[5])/((1-ss.tau)*(ss.phi_obar-par.phi_ubar)*ss.w)*100

        # c. savings
        moms['A_hh'] = ss.A_hh

        if vec:

            moms['A_hh_vec'] = np.zeros(par.Nfix)
            for i_fix in range(par.Nfix):

                if np.sum(ss.Dbeg[i_fix]) > 0.0:
                    moms['A_hh_vec'][i_fix] = np.sum(ss.a[i_fix]*ss.D[i_fix])/np.sum(ss.D[i_fix])
                else:
                    moms['A_hh_vec'][i_fix] = np.nan

        # d. print
        if do_print:
            for k,v in moms.items():

                if np.isscalar(v):                    
                    print(f'{k:25s} = {v:8.1f}',end='')        
                    k_vec = f'{k}_vec'
                    if vec and k_vec in moms:
                        print(f' | vec: [',end='')
                        for v_ in moms[k_vec]:
                            print(f'{v_:.1f} ',end='')
                        print(f']')
                    else:
                        print('') 

            print(f'moments in ss calculated in {elapsed(t0)}')     

    def calc_moms_path(self,do_print=False):
        """ calculate moments for calibration """

        par = self.par
        ss = self.ss
        path = self.path

        moms = self.moms

        # a. IRFs
        delta_IRF,_ = self.get_IRF('delta')
        lambda_u_IRF,_ = self.get_IRF('lambda_u')

        # b. moments
        delta = ss.delta + delta_IRF/100
        lambda_u = ss.lambda_u + lambda_u_IRF/100

        u_approx = 100*delta/(delta+lambda_u) - 100*ss.u
        u_approx_EU = 100*delta/(delta+ss.lambda_u) - 100*ss.u

        moms['w_share'] = 100*ss.w/(ss.px*ss.shock_TFP)
        moms['var_u'] = np.sum((100*(path.u[:,0]-ss.u))**2)

        if np.any(path.w[:,0] <= 0):
            moms['std_W'] = np.nan
        else:
            moms['std_W'] = np.sqrt(np.sum(((np.log(path.w[:,0])-np.log(ss.w)))**2))
        
        moms['timeshift'] = np.argmax(np.abs(lambda_u_IRF))-np.argmax(np.abs(delta_IRF))       
        moms['timeshift_obj'] = 100*(lambda_u_IRF[par.timeshift_target]-lambda_u_IRF.min())
        moms['EU_share'] = 100*np.cov(u_approx,u_approx_EU)[0,1]/np.var(u_approx)

        # d. print
        if do_print:
            for k,v in moms.items():
                if np.isscalar(v): 
                    print(f'{k:25s} = {v:12.4f}')        

    def calibrate_path_obj(self,x,skip_hh=False,check=True,do_print=False,do_print_full=False):
        """ evaluate objective function for calibration of path """

        t0 = time.time()

        par = self.par

        par.psi = x[0]
        par.xi = x[1]
        par.w_ss = x[2]
        
        psi = par.psi
        xi = par.xi
        w_ss = par.w_ss

        try:
            
            # a. find ss
            if skip_hh:
                steady_state.find_ss_SAM(self,do_print_full)
            else:
                self.find_ss(calib_beta=False,do_print=do_print_full)
            if do_print_full: print('')

            # b. compute jacs
            self.compute_jacs(skip_hh=skip_hh,skip_shocks=True,inputs_hh_all=['delta','lambda_u_s','RealR_ex_post','tau','hh_div'],do_print=do_print_full)               
            if do_print_full: print('')
            
            # c. find path
            self.find_transition_path(shocks=['shock_TFP'],do_end_check=False,do_print=do_print_full)
            if do_print_full: print('')
            
            # d. calculate moments
            self.calc_moms_path(do_print=do_print_full)
            if do_print_full: print('')

            # e. calculate error
            var_u = self.moms['var_u']
            EU_share = self.moms['EU_share']
            timeshift = self.moms['timeshift']
            timeshift_obj = self.moms['timeshift_obj']

            self.calibrate_path_objval = 0.0
            self.calibrate_path_objval += (timeshift-par.timeshift_target)**2            
            self.calibrate_path_objval += (100*timeshift_obj)**2            
            self.calibrate_path_objval += (100*(var_u-par.var_u_target))**2
            self.calibrate_path_objval += 10*(EU_share-par.EU_share)**2
        
            obj = self.calibrate_path_objval

            if do_print: 
                print(f'{psi = :6.4f}, {xi = :6.4f}, {w_ss = :6.4f}',end='')
                print(f' -> {EU_share = :5.2f}, {timeshift = :5.2f}, {timeshift_obj = :5.2f}, {var_u = :7.4f} -> {obj = :12.8f} [{elapsed(t0)}]')
            
            if check and self.calibrate_path_objval < par.eps_calib_path:
                raise StopIteration

        except StopIteration:

            raise

        except Exception as e:

            self.calibrate_path_objval = np.inf
            if do_print: print(f'{psi = :6.4f}, {xi = :6.4f}, {w_ss = :6.4f}: {e} [{elapsed(t0)}]')        

    def calibrate_path_ini(self,do_print=False,min_iter=500,max_iter=1000):
        """ find initial values for calibration of path """
        
        if do_print: print('find initial values:')
        t0 = time.time()

        x0 = np.nan*np.ones(3)
        bounds = ((0.1,3.0),(0.0001,0.1),(0.4,(self.par.epsilon_p-1)/self.par.epsilon_p-0.01))
        
        objval = np.inf
        found = False
        it = 0
        while True:
            
            if it%50 == 0:
                x0_ = np.array([np.random.uniform(bounds[i][0],bounds[i][1]) for i in range(3)])
                skip_hh = False
            else:
                x0_[:2] = np.array([np.random.uniform(bounds[i][0],bounds[i][1]) for i in range(2)])
                skip_hh = True

            self.calibrate_path_obj(x0_,skip_hh=skip_hh,check=False,do_print=do_print)
            objval_ = self.calibrate_path_objval

            if objval_ < objval:

                objval = objval_
                x0 = x0_
                found = True

            it += 1

            if found and it > min_iter: break
            if it > max_iter: raise ValueError('could not find initial value')

        if do_print: print(f'initial values found in {elapsed(t0)}')
        return x0

    def calibrate_path(self,x0,do_print=False,method='Nelder-Mead',maxiter=500):
        """ calibrate path """""

        par = self.par

        if do_print: print(f'calibrate path:')
        t0 = time.time()

        try:
            
            def obj_func(x): 
                self.calibrate_path_obj(x,do_print=do_print)
                return self.calibrate_path_objval

            bounds = ((0.1,10.0),(0.0001,0.5),(0.01,0.99))
            res = optimize.minimize(obj_func,x0,bounds=bounds,method=method,options={'maxiter':maxiter})        
            
            obj_func(res.x)
            objval = f'{self.calibrate_path_objval:12.8f}'

        except StopIteration:

            print('calibration stopped')
            objval = f'{self.calibrate_path_objval:12.8f} < {par.eps_calib_path:.1f}'

        if do_print:

            print('')
            print(f'{par.psi = :6.4f}')
            print(f'{par.xi = :6.4f}')
            print(f'{par.w_ss = :6.4f}')
            print(f'calibration obj.: {objval}')            
            print(f'done in {elapsed(t0)}')
    
    def calibrate_path_rep(self,x0,do_print=False,maxiter=500,min_iter_ini=50,max_iter_ini=500,max_iter_rep=5):
        """ calibrate path, repeat until convergence"""
        
        par = self.par

        it = 0
        while True:
            
            self.calibrate_path(x0,do_print=do_print,maxiter=maxiter,min_iter_ini=min_iter_ini,max_iter_ini=max_iter_ini)

            if self.calibrate_path_objval < par.eps_calib_path: 

                break

            else:

                # not converged
                objval = self.calibrate_path_objval
                print('')
                print(f'{objval = :.1f}, not converged, repeat')
                print('')

                # new initial guess
                x0 = self.calibrate_path_ini()

            it += 1
            if it > max_iter_rep: break
    
    ##########
    # policy #
    ##########

    def find_policy(self,pol,du_target=None,basepol_value=None,basepol_rho=None):
        """ find policy - base or unemployment target """

        model = self.copy()

        # a. set shock or target
        if du_target is None:
        
            assert np.isscalar(basepol_value)
            assert np.isscalar(basepol_rho)
            assert du_target is None

            K = 120
            shocks = {f'd{pol}':np.zeros(model.par.T)}
            shocks[f'd{pol}'][:K] += basepol_value*basepol_rho**np.arange(K)    

            model.H_U = self.H_U # copy can remove it due to set_macro_auto

        else:
            
            assert basepol_value is None
            assert basepol_rho is None
            assert not du_target is None

            model.shocks = [shock for shock in model.shocks if not shock == pol]
            model.unknowns += [pol]
            model.targets += ['errors_u_target']

            model.allocate_GE(update_hh=False,ss_nan=False)
            model.compute_jacs(skip_hh=True,skip_shocks=True)

            shocks = {'du_target':du_target}
        
        # b. find policy
        model.find_transition_path(shocks=shocks,do_end_check=False)
        model.moms = {}
        model.calc_moms_path()
        
        return model