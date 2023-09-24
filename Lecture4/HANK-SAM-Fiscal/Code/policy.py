from IPython.display import display

from types import SimpleNamespace
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from figures import create_fig, save_fig

def create_models(basepol,model_specs,do_print=True):
    """ create models and run calibrations """

    for name,specs in model_specs.items():

        if do_print: print(f'### {name} ###')

        model = basepol.copy(name=name)

        # a. set parameters
        for k,v in specs['pardict'].items():
            model.par.__dict__[k] = v
        
        # b. update macros and grids
        model.set_macros_auto()
        model.create_beta_grid()

        # c. save
        model.save() # might be overwritten later
        model.compress()

        # d. calibrate
        if 'calib_beta' in model_specs[name]:
            calib_beta = model_specs[name]['calib_beta']
        else:
            calib_beta = False
                
        if 'recalib_beta' in model_specs[name]:
            recalib_beta = model_specs[name]['recalib_beta']
        else:
            recalib_beta = False

        if 'calib_path' in model_specs[name]:
            calib_path = model_specs[name]['calib_path']
        else:
            calib_path = False

        if 'calib_path_ini' in model_specs[name]:
            calib_path_ini = model_specs[name]['calib_path_ini']
        else:
            calib_path_ini = False            

        if any([calib_beta,recalib_beta,calib_path]):
            
            if do_print: print('')
            try:
                model.full_run(calib_beta=calib_beta,recalib_beta=recalib_beta,calib_path=calib_path,calib_path_ini=calib_path_ini,
                                do_print_full=True,do_print=True,do_save=True)        
            except Exception as e:            
                print(f'failed: {e}')

def add_var_u_G(model,do_print=True):
    """ var_u for G shock """

    if do_print: print(f' [var_u_G] ',end='')

    try:
        
        # a. compute
        model_demand = model.copy()
        model_demand.H_U = model.H_U # copy can remove it due to set_macro_auto
        model_demand.find_transition_path(shocks=['G'],do_end_check=False)
        model_demand.calc_moms_path()

        # b. save
        model.moms['var_u_G'] = model_demand.moms['var_u']

        if do_print: print('')

    except Exception as e:

        if do_print: print(f' failed, {e}')               

def find_policy_and_save(model,pol,du_target,do_print=True):
    """ find policy and save """

    try:
        polmodel = model.find_policy(pol,du_target=du_target)
        polmodel.name = model.name + f'_{pol}'
        save_results(polmodel)
        if do_print: print('')
    except Exception as e:
        if do_print: print(f' failed, {e}')

def save_results(model):
    """ save results """

    results = {}
    results['par'] = model.par.__dict__
    results['ss'] = {k:model.ss.__dict__[k] for k in model.varlist}
    results['path'] = {k:model.path.__dict__[k] for k in model.varlist}
    results['moms'] = model.moms

    with open(f'saved/{model.name}.p', 'wb') as f:
        pickle.dump(results,f)

def policies(model,pols,model_flex=None,model_RA=None,pols_flex=None,pols_RA=None,
            basepol='G',basepol_value_rel_ss=0.001,basepol_rho=0.965,do_print=True):
    """ find policies """
    
    do_flex = model_flex is not None
    do_RA = model_RA is not None

    assert basepol in pols
    if do_flex: assert not pols_flex is None
    if do_RA: assert not pols_RA is None

    # a. baseline policy
    if do_print: print(f'{basepol = }',end='')
    try:
        basepol_value = basepol_value_rel_ss*model.ss.__dict__[basepol]
        basepolmodel = model.find_policy(basepol,basepol_value=basepol_value,basepol_rho=basepol_rho)
        basepolmodel.name = model.name + f'_{basepol}'
        save_results(basepolmodel)
        if do_print: print('')
    except Exception as e:
        if do_print: print(f' failed, {e}\n')
        return
    
    # e. implied policy target
    du_target = basepolmodel.path.u[:,0]-basepolmodel.ss.u

    # e. other policies
    for pol in pols:

        # i. sticky prices
        if pol == basepol: 
            save_results(basepolmodel)
        else:
            if do_print: print(f' {pol}',end='')
            find_policy_and_save(model,pol,du_target,do_print=True)

        # ii. flexible price
        if do_flex and pol in pols_flex:
            if do_print: print(f' {pol}_flex',end='')
            find_policy_and_save(model_flex,pol,du_target,do_print=True)

        # iii. RA
        if do_RA and pol in pols_RA:
            if do_print: print(f' {pol}_RA',end='')
            find_policy_and_save(model_RA,pol,du_target,do_print=True)

def run_all(basemodel,model_names,pols,pols_flex=None,pols_RA=None,
            basepol='G',basepol_value_rel_ss=0.001,basepol_rho=0.965,do_print=True):
    
    for name in model_names:

        if do_print: print(f'### {name} ###\n')

        # a. load
        model = basemodel.copy(name=name)
        model.name = name
        model.load()

        # b. full run
        print('full run',end='')
        try:
            model.full_run()
            add_var_u_G(model,do_print=do_print)  
            save_results(model)
        except Exception as e:
            print(f' failed, {e}')           

        # c. other models
        do_flex = not pols_flex is None
        if do_flex:

            print(f'_flex')
            try:
                model_flex = model.copy()
                model_flex.name = model.name + '_flex'
                model_flex.par.phi = 0.001
                model_flex.full_run(do_ss=False,skip_hh=True)
                save_results(model_flex)      
            except Exception as e:
                print(f' failed, {e}') 
                do_flex = False
        else:
            model_flex = None

        do_RA = (not pols_RA is None) and np.isclose(model.par.div_hh,1.0)
        if do_RA:
            
            print(f'_RA',end='')
            try:
                model_RA = model.copy()
                model_RA.name = model.name + '_RA'
                model_RA.par.RA = True
                model_RA.full_run(do_ss=False,skip_hh=True)
                add_var_u_G(model_RA)
                save_results(model_RA)
            except Exception as e:
                print(f' failed, {e}') 
                do_RA = False    

        else:
            model_RA = None            

        # d. baseline policy
        policies(model,pols,model_flex=model_flex,model_RA=model_RA,pols_flex=pols_flex,pols_RA=pols_RA,
                basepol=basepol,basepol_value_rel_ss=basepol_value_rel_ss,basepol_rho=basepol_rho,do_print=do_print)

        if do_print: print('')

def load_model(name):
    """ load model """

    with open(f'saved/{name}.p', 'rb') as f:
        results = pickle.load(f) 

    model = SimpleNamespace()    
    model.par = SimpleNamespace(**results['par'])
    model.ss = SimpleNamespace(**results['ss'])
    model.path = SimpleNamespace(**results['path'])
    model.moms = results['moms']

    return model

def load_models(model_names,flex=False,RA=False,do_print=True):
    """ load models """

    model_names_all = [name for name in model_names]
    if flex: model_names_all += [f'{name}_flex' for name in model_names]
    if RA: model_names_all += [f'{name}_RA' for name in model_names]

    models = {}
    for name in model_names_all:
        try:
            models[name] = load_model(name)
        except:
            pass
    
    return models

def load_pol_models(model_names,pols,pols_flex=None,pols_RA=None,do_print=True):
    """ load policy models """

    model_names_all = [name for name in model_names]
    pols_all = [pol for pol in pols]

    flex = not pols_flex is None
    if flex: model_names_all += [f'{name}_flex' for name in model_names]

    RA = not pols_RA is None
    if RA: model_names_all += [f'{name}_RA' for name in model_names]

    pol_models = {}
    for name in model_names_all:
        for pol in pols_all:

            fullname = f'{name}_{pol}'

            try:
                pol_models[(name,pol)] = load_model(fullname)
            except:
                pass
    
    return pol_models

def pol_figures(name,pols,varlist,pol_models):
    """ plot policy figures """

    # a. label
    var_labels = {
        'G':'government spending',
        'u_bar':'unemployment duration',
        'u':'unemployment rate',
        'tau':'tax rate',   
        'qB':'value of government debt',
    }
      
    pol_labels = {
        'G':'public spending',
        'hh_transfer':'household transfer',
        'retention_subsidy':'retention subsidy',
        'hiring_subsidy':'hiring subsidy',
        'phi_obar':'UI level',
        'u_bar':'UI duration',
        }

    pol_ls = {
        'G':'-',
        'hh_transfer':'-',
        'retention_subsidy':'--',
        'hiring_subsidy':'--',
        'phi_obar':':',
        'u_bar':':',
        }
    
    pol_lw = {
        'G':'1',
        'hh_transfer':'1',
        'retention_subsidy':'1',
        'hiring_subsidy':'1',
        'phi_obar':'2',
        'u_bar':'2',
        }    

    # b. figure for each variable
    varlist_pct = ['G','u','tau']
    for varname in varlist:

        fig,ax = create_fig(figsize=(8,5))

        for pol in pols:

            if (name,pol) in pol_models:

                model = pol_models[(name,pol)]
                
                ls = pol_ls[pol]
                lw = pol_lw[pol]
                pol_label = pol_labels[pol]

                if varname in varlist_pct:
                    ax.plot((model.path.__dict__[varname][:,0]/model.ss.__dict__[varname]-1)*100,label=pol_label,ls=ls,lw=lw)
                else:
                    ax.plot(model.path.__dict__[varname][:,0]-model.ss.__dict__[varname],label=pol_label,ls=ls,lw=lw)

        if varname in varlist_pct:
            ylabel = '%'
        else:
            ylabel = ''

        save_fig(fig,ax,filename=f'pol_{name}_{varname}',
            title=var_labels[varname],ylabel=ylabel,legend=varname in ['u'],ncol=3,T_max=48)
        
        plt.show()
    
def moms(models,momnames=None,model_labels=None,tex=None):
    """ table of moments """

    if momnames is None: momnames = ['MPC_qtr','C_drop_ss','C_drop_ex','EU_share','timeshift','A_hh','var_u']

    # a. prep    
    Nmoms = len(momnames)
    Nmodels = len(models)

    # b. moments
    M = np.nan*np.ones((Nmoms,Nmodels))
    for j,(modelname,model) in enumerate(models.items()):
        for i,momname in enumerate(momnames):
            if momname in model.moms:
                M[i,j] = model.moms[momname]

    # c. table
    mom_labels = {
        'MPC_qtr':'MPC (qtr.)',
        'C_drop_ss':'C. drop in unemployment',
        'C_drop_ex':'Relative C. drop at UI exit',
        'EU_share':'Share of EU',
        'timeshift':'Time-shift',
        'var_u':'TFP',
        'var_u_G':'G',
        'A_hh':'A_hh', 
    }  

    df = pd.DataFrame()

    if model_labels is None:
        df[''] = [modelname for modelname in models.keys()]
    else:
        df[''] = [model_labels[modelname] for modelname in models.keys()]

    for i,momname in enumerate(momnames): 
        df[mom_labels[momname]] = ['']*Nmodels
        for j in range(Nmodels):
            if np.isnan(M[i,j]): continue
            if 'var_u' in momname:
                df[mom_labels[momname]].values[j] = f'{M[i,j]:.4f}'
            else:
                df[mom_labels[momname]].values[j] = f'{M[i,j]:.2f}'

    display(df)

    if not tex is None:
        df = df.rename(columns={model_labels[modelname]:f'\\textbf{{{model_labels[modelname]}}}' for modelname in models.keys()})
        df.style.hide(axis="index").to_latex(f'results/{tex}.tex',hrules=True)

    return df

def fiscal_multipliers(pol_models,
                       select=None,
                       tex=None,
                       do_rel=False,not_G_rel=False,model_labels=None,nan_value=''):
    """ table of fiscal multipliers """

    if not select is None:
        pol_models = {k:v for k,v in pol_models.items() if k[0] in select}

    # a. prep    
    names = []
    pols = []
    for k in pol_models.keys():
        if not k[0] in names: names.append(k[0])
        if not k[1] in pols: pols.append(k[1])

    Nnames = len(names)
    Npols = len(pols)

    # b. multiplier
    M = np.nan*np.ones((Npols,Nnames))
    for i,pol in enumerate(pols):
        for j,name in enumerate(names):

            if (name,pol) in pol_models:
                
                model = pol_models[(name,pol)] 
                tax_expenditures = np.sum(model.path.__dict__['tau'][:,0]*model.path.__dict__['Yt_hh'][:,0]-model.ss.__dict__['tau']*model.ss.__dict__['Yt_hh'])
                unemployment = np.sum(model.path.__dict__['u'][:,0]-model.ss.__dict__['u'])

                M[i,j] = -unemployment/tax_expenditures

    if do_rel:
        for i,pol in enumerate(pols):
            if pol == 'G':
                if not_G_rel:
                    orig = M[i,:].copy()
                    M /= M[i,:]
                    M[i,:] = orig
                else:
                    M /= M[i,:]
            
    # c. table
    df = pd.DataFrame()

    pol_labels = {
        'G':'G',
        'hh_transfer':'transfer',
        'retention_subsidy':'retention',
        'hiring_subsidy':'hiring',
        'phi_obar':'UI level',
        'u_bar':'UI duration',
    }  

    if model_labels is None:
        df[''] = names
    else:
        df[''] = [model_labels[name] for name in names]

    for i,pol in enumerate(pols):
        df[pol_labels[pol]] = ['']*len(names)
        for j in range(len(names)):
            if np.isnan(M[i,j]):
                df[pol_labels[pol]].values[j] = nan_value
            else:
                if do_rel and not_G_rel and pol == 'G':
                    df[pol_labels[pol]].values[j] = f'1.0 [{M[i,j]:.2f}]'
                else:
                    df[pol_labels[pol]].values[j] = f'{M[i,j]:.2f}'
                    
    display(df)

    if not tex is None:
        df = df.rename(columns={pol_labels[pol]:f'\\textbf{{{pol_labels[pol]}}}' for pol in pols})
        df.style.hide(axis="index").to_latex(f'results/fiscal_multipliers_{tex}.tex',hrules=True)