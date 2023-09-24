
import matplotlib.pyplot as plt
import numpy as np
from IHANKModel import HANKModelClass
from scipy.interpolate import CubicSpline
from copy import deepcopy
from GEModelTools import lag, simulate_hh_path 
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import norm

import matplotlib.pyplot as plt   
from seaborn import color_palette, set_palette
from matplotlib import rc

plt.style.use('seaborn-white')
set_palette("colorblind")

rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{underscore}')

############
# settings #
############

abs_value = []

pctp = ['iF_s','piF_s','piF','pi','piNT','piH','ppi','r','i',
        'r_ann','pi_ann','Q','Dompi','di','NFA', 'r_NFA', 'piH_ann']

paths_defaults = {
    'standard':
        ['Y','C','I','Exports','pi_ann',
         'YNT','CNT','N','Imports','r_ann',
         'YT','CT','labor_comp','ToT','Q'],
    'standard_vs_data':
        ['Y','C','I','Exports','pi_ann',
         'YNT','CNT','N','Imports','r_ann',
         'YT','CT','labor_comp','piH_ann','Q']
}
pathlabels = {

    'Y_s':'Foreign output ($Y^*$)',
    'C_s':'Foreign demand ($C^*$)',
    'piF_s':'Foreign inflation ($\pi^*_F$)',
    'iF_s':'Foreign interest rate ($i^*$)',
    'Y':'GDP ($Y$)',
    'YNT':'Non-Tradeable VA ($Y_{NT}$)',
    'YT':'Tradeable VA ($Y_{T}$)',
    'C':'Consumption ($C$)',
    'CNT':'Consumption of non-tradeables ($C_{NT}$)',
    'CT':'Consumption of tradeables ($C_{T}$)',
    'I':'Investment ($I$)',
    'N':'Employment ($N$)',
    'labor_comp':'Wage Income ($wN$)',
    'Q':'Real exchange rate ($Q$)',
    'ppi':'PPI ($\pi^{PP}$)',
    'pi':'Inflation ($\pi$)',
    'pi_ann':'Inflation ($\pi$) (annual)',
    'r':'Real interest rate ($r$)',
    'r_ann':'Real interest rate ($r$) (annual)',
    'Exports':'Exports',
    'Imports':'Imports',
    'NX':'Net exports',
    'piH_ann':'$P_H$ Inflation ($\pi_H$) (annual)',
    'ToT':'Terms of Trade',

    'tau':'Tax rate ($\tau$)',
    'B':'Government debt ($B$)',

    'eps_beta':'Discount factor ($\\beta$)',

    'PT_PNT':r'$P_{T}/P_{NT}$',
    
    'rel_C_hh':r'Sectoral consumption, $C_T^{hh} / C_{NT}^{hh}$',
    'rel_wn_hh':r'Sectoral income, $w_T N_T / (w_{NT}N_{NT}$)',

}


################
# steady state #
################

def plot_MPCs(model,model_alt=None,label='HANK',label_alt='Two-Asset HANK',plot_RANK=True):
    """ Plot MPCs """

    # MPCs in data (Fagereng et al. 2021: "MPC heterogeneity and household balance sheets")
    MPC_data = {}
    MPC_data['x'] = np.arange(6)
    MPC_data['y'] =  model.par.MPC_est
    
    cs = CubicSpline(MPC_data['x'],MPC_data['y'] )
    Nquarters = 11 
    
    MPC_data['x_int'] = np.arange(Nquarters)/4
    MPC_data['x_int_Q'] = np.arange(Nquarters)
    
    y_int = cs(MPC_data['x_int'])/4
    MPC_data['y_int'] = y_int / np.sum(y_int[:4]) * MPC_data['y'][0] 

    # Model MPCs
    path_org = deepcopy(model.path)
    ann_MPC,MPCs = model.ann_MPCs()
    model.path = path_org

    if model_alt is not None:

        path_org = deepcopy(model_alt.path)
        ann_MPC_alt,MPCs = model_alt.ann_MPCs()
        model_alt.path = path_org      

    # RANK
    beta = 1/(1+model.ss.ra)
    MPC_RANK = (1-(beta*(1+model.ss.r)**(1/model.par.CRRA))/(1+model.ss.r))*beta**np.arange(model.par.T)
    ann_MPC_RANK = np.zeros(round(model.par.T/4))
    for j in range(round(model.par.T/4)):
        ann_MPC_RANK[j] = np.sum(MPC_RANK[j*4:(1+j)*4])   

    # Plot
    lsize = 1.8
    sizetupple = (4.4*1,2.9*1)
    
    fig = plt.figure(figsize=sizetupple)
    ax = fig.add_subplot(1,1,1)

    ax.plot(np.zeros(6), '-', color='black')
    ax.plot((ann_MPC[:6]), '-', label=label, color='C2', linewidth=lsize)
    if model_alt is not None: ax.plot((ann_MPC_alt[:6]), '-.', label=label_alt,color='C3',linewidth=lsize)
    if plot_RANK: ax.plot(ann_MPC_RANK[:6], '--', label='RANK', color='C1', linewidth=lsize)
    ax.plot(MPC_data['y'], linestyle='None', marker='o', label='Fagereng et al. (2021)', color='C0')
    
    ax.set_xlabel('Years', fontsize=16)
    ax.set_ylabel('Annual MPC', fontsize=12)
    ax.set_xlim([0, 6-1])
    
    ax.legend(frameon=True, fontsize=12)
    
    fig.tight_layout()    
    
    return fig 

def PE_MonPol_shock(model,do_KMV=False):
    """ Plot consumption response to interest rate change"""
    
    # Kaplan et al. (2018) shock -  exact sources?
    x_KMV_shock = np.array([0.02, 0.3,0.8,1.2,2.0,2.4,3.2,3.9,4.7,5.1,6.1,6.6,7.6,8.5,9.2,10,11,12,13])
    y_KMV_shock = np.array([-0.666,-0.582,-0.453,-0.360,-0.252,-0.214,-0.145,-0.107,-0.072,-0.055,-0.035,-0.036,-0.024,-0.019,-0.020,0,0,0,0])/4 # convert to quarterly 
    
    cs_shock = CubicSpline(x_KMV_shock, y_KMV_shock)
    y_shock_int = cs_shock(np.arange(0,10))
    
    KMV_shock = np.zeros(model.par.T)
    KMV_shock[:10] = y_shock_int

    x_KMV =   np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
    y_KMV = np.array([0.14*1000,62,19,3,-8,-14,-18,-19,-19,-18,-16,-13,-13,-13,-12,-12,-12,-12,-11,-11,-11])/1000
    
    cs = CubicSpline(x_KMV, y_KMV)
    y_int = cs(np.arange(0,21))
     
    dra = KMV_shock / 100 
    
    # RANK
    beta = 1/(1+model.ss.ra)
    C_RANK = np.zeros(model.par.T)        
    for tt in range(model.par.T):
        t = model.par.T-1-tt
        if t == model.par.T-1:
            C_RANK[t] = model.ss.C            
        else:
            C_RANK[t] =  model.ss.C-1/model.par.CRRA * beta * np.sum(dra[t+1:]) + model.ss.A*(1-beta)*dra[0]
            
    # Model
    C_ss = np.sum(model.ss.c*model.ss.D)
    model._set_inputs_hh_all_ss()

    model.path.ra[:,0] = model.ss.ra + dra
    model.solve_hh_path()
    model.simulate_hh_path()    
        
    dC = np.zeros(model.par.T)*np.nan
    for t in range(model.par.T):
        dC[t] = deepcopy(np.sum((model.path.c[t] * model.path.D[t])) - C_ss) 

    model._set_inputs_hh_all_ss()

    # plot
    t = 21
    xaxis = np.arange(0,t)
    lsize = 1.8 
    set_palette("colorblind")
    
    fig = plt.figure(figsize=(4.7*1,3.2*1))
    ax = fig.add_subplot(1,1,1)

    ax.plot(xaxis,np.zeros(t), '--', color='black')
    ax.plot(xaxis,(dC[:t]/model.ss.C)*100, label='HANK', color='C2', linewidth=lsize)
    ax.plot(xaxis,(C_RANK[:t]/model.ss.C-1)*100, '-', label='RANK', color='C1', linewidth=lsize)
    if do_KMV: ax.plot(xaxis,y_int, '--', label='Kaplan et al. (2018)', color='C0', linewidth=lsize)
    
    ax.set_xticks(np.arange(0, t, 4))
    ax.set_xlabel('Quarters', fontsize=16)
    ax.set_ylabel('\% diff. to s.s.', fontsize=12)
    ax.set_xlim([0, t-1])
    
    ax.legend(frameon=True, prop={'size': 12})
    
    fig.tight_layout()        
    
    return fig

def plot_jac_columns(model_RA,model_HA,model_HA2A=None,calc_jac=True):
    """ Plot Jacobian columns """
    
    # Compute Jacobians
    if calc_jac:
        model_HA._compute_jac_hh()
        if model_HA2A is not None:
            model_HA2A._compute_jac_hh()
            HA2A_dC_dI = model_HA2A.jac_hh[('C_hh', 'UniformT')]
            HA2A_dC_dr = model_HA2A.jac_hh[('C_hh', 'ra')]   

    HA_dC_dI = model_HA.jac_hh[('C_hh', 'UniformT')]
    HA_dC_dr = model_HA.jac_hh[('C_hh', 'ra')]

    # Convert RA_dC_dr to ra dating 
    if not model_RA is None:
        
        RA_dC_dI = model_RA.par.M_Y.copy() 
        RA_dC_dr = model_RA.par.M_R.copy() 
    
        RA_dC_dr[:,1:] = RA_dC_dr[:,:-1].copy()
        RA_dC_dr[:,0] = RA_dC_dI[:,0] * model_RA.ss.r

    # Plot
    set_palette("colorblind")
    lsize = 1.5
    Nquarters = 30 

    columns = [0,4, 8, 12, 16, 20, 24]
    alphalist = np.flip(np.linspace(0.4,1,len(columns)))

    l1,l2,l3 = 'HANK', 'RANK', 'Two-asset HANK'

    fig = plt.figure(figsize=(8,3.0))

    ax = fig.add_subplot(1,2,1)
    ax.set_title(r'$\mathbf{M}$')

    l1_c,l2_c,l3_c = l1, l2, l3
    for col in columns:   
        ax.plot(HA_dC_dI[:Nquarters, col],'-',label=l1_c,color='C2',linewidth=lsize)
        if model_HA2A is not None:
            ax.plot(HA2A_dC_dI[:Nquarters, col], '-.', label=l3_c, color='C3', linewidth=lsize)
        
        if not model_RA is None:
            ax.plot(RA_dC_dI[:Nquarters, col], '--', label=l2_c, color='C1', linewidth=lsize)    

        l1_c, l2_c, l3_c = '_nolegend_', '_nolegend_', '_nolegend_' # suppress extra legends

    ax.plot(np.zeros(Nquarters), '-', color='black')

    ax.set_xticks(np.arange(0, Nquarters, 4))
    ax.set_xlim([0, Nquarters-1])
    ax.set_xlabel('Quarters', fontsize=16)
    ax.set_ylabel(f'$dC$', fontsize=12)

    ax = fig.add_subplot(1,2,2)
    ax.set_title(r'$\mathbf{M^r}$')

    l1_c,l2_c,l3_c = l1, l2, l3
    for i,col in enumerate(columns):   
        ax.plot(HA_dC_dr[:Nquarters, col],'-',label=l1_c,color='C2',linewidth=lsize, alpha=alphalist[i])
        if model_HA2A is not None:
            ax.plot(HA2A_dC_dr[:Nquarters, col], '-.', label=l3_c, color='C3', linewidth=lsize, alpha=alphalist[i])
        if not model_RA is None:
            ax.plot(RA_dC_dr[:Nquarters, col], '--', label=l2_c, color='C1', linewidth=lsize, alpha=alphalist[i])    
        l1_c, l2_c, l3_c = '_nolegend_', '_nolegend_', '_nolegend_' # suppress extra legends

    ax.plot(np.zeros(Nquarters), '--', color='black')
    
    ax.set_xticks(np.arange(0, Nquarters, 4))
    ax.set_xlim([0, Nquarters-1])
    ax.set_xlabel('Quarters', fontsize=16)
    ax.legend(frameon=True, fontsize=12)
    ax.set_ylabel(f'$dC$', fontsize=12)

    plt.tight_layout()

    return fig 

########
# IRFs #
########

def get_dX(varname, model, absvalue=False, scaleval=1.0):

    pathvalue = getattr(model.path,varname)  
    ssvalue = getattr(model.ss,varname)

    if absvalue:
        dX = (pathvalue-ssvalue) * scaleval
    else:
        dX = (pathvalue-ssvalue) * scaleval / ssvalue

    return dX

def get_single_HA_IRF(model,inputvar,outpuvar,scaleval):

    ss_outputvar = getattr(model.ss,outpuvar)

    dX = getattr(model.path,inputvar)
    X_ss = getattr(model.ss, inputvar)

    X_jac_hh = model.jac_hh[(outpuvar,inputvar)]
    
    return X_jac_hh @ (dX-X_ss)*scaleval / ss_outputvar * 100 

def C_decomp_HA_v_RA(modellist,T_max,lwidth,disp_income=True,scale=True,test_plot=False):   
    """ Decomposition of consumption IRFs for HA and RA-IM models """

    C_decomp = {}

    for i,model in enumerate(modellist):
        hh = model.par.HH_type

        if scale:
            scaleval = getattr(model.par,'scale')
        else:
            scaleval = 1 

        if model.par.HH_type in ['HA', 'HA-2A']:

            ss_C = getattr(model.ss, 'C')
            
            C_decomp[hh] = {}
            C_decomp[hh]['Total'] = get_dX('C',model,scaleval=scaleval)*100
            C_decomp[hh]['Labor income'] = get_single_HA_IRF(model,'wnNT','C_hh',scaleval) + get_single_HA_IRF(model,'wnT','C_hh',scaleval)
            C_decomp[hh]['Taxes'] = get_single_HA_IRF(model,'tau','C_hh',scaleval)
            
            X_ss = getattr(model.ss,'ra')
            dr = getattr(model.path,'ra') - X_ss
            dr[0] = 0.  
            ra_jac_hh = model.jac_hh[('C_hh','ra')]
            C_decomp[hh]['r'] = ra_jac_hh @ dr*scaleval / ss_C * 100   
            
            dra   = getattr(model.path, 'ra') - X_ss
            dra[1:] = 0 
            C_decomp[hh]['ra'] = ra_jac_hh @ dra*scaleval / ss_C * 100   

            if 'rai' in model.inputs_hh:
                drai   = getattr(model.path, 'rai') - model.ss.rai 
                drai[1:] = 0.
                rai_jac_hh = model.jac_hh[('C_hh','rai')]
                C_decomp[hh]['ra'] += rai_jac_hh @ drai*scaleval / ss_C * 100   
                C_decomp[hh]['r'] += rai_jac_hh @ dr*scaleval / ss_C * 100   

            C_decomp[hh]['test'] = C_decomp[hh]['ra'] + C_decomp[hh]['r']  \
                                    + C_decomp[hh]['Labor income'] +  C_decomp[hh]['Taxes'] 

        elif model.par.HH_type == 'RA-IM':

            # RA 
            ss_C = getattr(model.ss, 'C')

            C_decomp[hh] = {}
            C_decomp[hh]['Total'] = get_dX('C',model,scaleval=scaleval) *100

            sT = model.par.sT
            dI = get_dX('wnNT',model,scaleval=scaleval,absvalue=True)*(1-sT) + get_dX('wnT',model,scaleval=scaleval,absvalue=True)*sT
            C_decomp[hh]['Labor income'] = model.par.M_Y @  dI / ss_C * 100 
            C_decomp[hh]['Taxes'] = model.par.M_Y @ get_dX('tau', model,scaleval=scaleval,absvalue=True) / ss_C * 100 

            X_ss = getattr(model.ss,'r')
            dX = getattr(model.path,'r')  - X_ss
            C_decomp[hh]['r'] = model.par.M_R @ dX*scaleval / ss_C * 100   

            dra = getattr(model.path,'ra')- X_ss
            dra[1:] = 0 
            C_decomp[hh]['ra'] = model.par.M_Y @ dra*scaleval / ss_C * 100   

            C_decomp[hh]['test'] = C_decomp[hh]['ra'] + C_decomp[hh]['r']  \
                                + C_decomp[hh]['Labor income'] +  C_decomp[hh]['Taxes'] 
        else:
            raise ValueError('HH type need to be either HA or RA-IM')

    # Plot
    ncols,nrows = 2,1
    set_palette("colorblind")

    fig = plt.figure(figsize=(4.3*ncols,3.6*nrows))
    
    for i,model in enumerate(modellist):
        hh = model.par.HH_type 

        ax = fig.add_subplot(nrows,ncols,i+1)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_xticks(np.arange(0, T_max, 4))

        ax.plot(np.zeros(T_max), '-', color='black')   
        
        if disp_income:
            dY = C_decomp[hh]['Labor income'][:T_max] + C_decomp[hh]['Taxes'][:T_max]
            ax.plot(np.arange(T_max),dY,label='Disp. income', linewidth=lwidth, linestyle='--')
        else:
            ax.plot(np.arange(T_max),C_decomp[hh]['Labor income'][:T_max],label='Labor income', linewidth=lwidth, linestyle='--')
            ax.plot(np.arange(T_max),C_decomp[hh]['Taxes'][:T_max],label='Taxes', linewidth=lwidth, linestyle=':')
        
        ax.plot(np.arange(T_max),C_decomp[hh]['r'][:T_max], label='Interest rate', linewidth=lwidth, linestyle='-.', color='C2')
        ax.plot(np.arange(T_max),C_decomp[hh]['ra'][:T_max], label='Revaluation', linewidth=lwidth, color='Firebrick', marker=',')     
        ax.plot(np.arange(T_max),C_decomp[hh]['Total'][:T_max], label='Total', linewidth=lwidth, marker='o', color='C3')
        
        if test_plot:
            ax.plot(np.arange(T_max),C_decomp[hh]['test'][:T_max], label='Test', linewidth=lwidth, linestyle='-.', color='black')

        ax.set_xlim([0, T_max-1])       
        ax.set_ylabel('$\%$ diff. to s.s.')
        ax.set_xlabel('Quarters')

        # if hh == 'RA':
        #     ax.set_title('RANK')
        # else:
        #     ax.set_title('HANK')
        ax.set_title(model.name)

    ax.legend(loc="best", frameon=True)  
    fig.tight_layout()
   
    return fig   

def trad_nontrad_decomp(models,labels=['HANK','RANK'],T_max=21):
    """ Plot decomposition of consumption import """

    if type(models) is dict: models = [model for model in models.values()]
    
    # Get IRFs
    resp_dict = {}    
    for model_ in models:
        
        name = model_.name
        scale = model_.par.scale
        resp_dict[name] = {}
        resp_dict[name]['T'] = {}
        resp_dict[name]['NT'] = {}

        etaT = model_.par.etaT

        resp_dict[name]['dC']   = (model_.path.C/model_.ss.C-1)*scale*100
        resp_dict[name]['T']['dC']   = (model_.path.CT/model_.ss.CT-1)*scale*100
        resp_dict[name]['NT']['dC']   = (model_.path.CNT/model_.ss.CNT-1)*scale*100
        resp_dict[name]['T']['dP'] = - etaT*(model_.path.PT/model_.path.P -1)*scale*100
        resp_dict[name]['NT']['dP'] = - etaT*(model_.path.PNT/model_.path.P -1)*scale*100
        
    parvalslabels = ['Non-tradeable C', 'Tradeable C']
        
    # Plot
    lwidth = 2.3
    ncols,nrows = 2,2
    fig = plt.figure(figsize=(4.3*ncols/1.1,3.6*nrows/1.2))

    temp = 0 
    for j,sec in enumerate(['NT', 'T']):    
        for i,(model_,label) in enumerate(zip(models,labels)):

            name = model_.name

            ax = fig.add_subplot(nrows,ncols,temp+1)
            ax.plot(np.zeros(T_max), '-', color='black')   

            ax.plot(np.arange(T_max),resp_dict[name][sec]['dC'][:T_max],label='Total', linewidth=lwidth, linestyle='-')
            ax.plot(np.arange(T_max),resp_dict[name]['dC'][:T_max],label='Scale effect', linewidth=lwidth, linestyle='--')
            ax.plot(np.arange(T_max),resp_dict[name][sec]['dP'][:T_max],label='Substitution effect', linewidth=lwidth, linestyle=':')
                  
            ax.set_xlim([0, T_max-1])
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.set_xticks(np.arange(0, T_max, 4))

            ax.set_ylabel('$\%$ diff. to s.s.')
            ax.set_xlabel('Quarters')
            
            ax.set_title(label)

            if temp==0: ax.legend(loc="best", frameon=True)  
            temp +=1

            if i == 0:
                pad=0.8
                ax.annotate(parvalslabels[j], xy=(-0.5, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size=18, ha='right', va='center', rotation=90)
    
    fig.tight_layout()
    
    return fig 

def _plot_IRFs(ax,model,pathname,scale,lstyle,color,lwidth,label,T_max):

    global abs_value, pctp, pathlabels

    if pathname in pathlabels:
        pathlabel = pathlabels[pathname]
    else:
        pathlabel = pathname

    if scale:
        scaleval = getattr(model.par,'scale')
    else:
        scaleval = 1 

    # ssvalue and pathvalue
    ssvalue = getattr(model.ss,pathname)  
    pathvalue = getattr(model.path,pathname)
    dpathvalue = pathvalue - ssvalue

    T_max = np.fmin(pathvalue.size,T_max)
    
    # plot
    if pathname in abs_value:                     
    
        ax.plot(np.arange(T_max),(dpathvalue[:T_max])*scaleval,label=label,linestyle=lstyle,color=color,linewidth=lwidth)
        ax.set_ylabel('abs. diff. to s.s.')
        ax.set_title(pathlabel)
    
    elif pathname in pctp:
    
        ax.plot(np.arange(T_max),100*(dpathvalue[:T_max])*scaleval,label=label,linestyle=lstyle,color=color,linewidth=lwidth)
        ax.set_ylabel('\%-points diff. to s.s.')
        ax.set_title(pathlabel)
    
    elif pathname == 'NX':   
        
        pathvalue_IM = getattr(model.path,'Imports')   
        pathvalue_EX = getattr(model.path,'Exports')  
        ssvalue_IM = getattr(model.ss,'Imports')
        ssvalue_EX = getattr(model.ss,'Exports')
        dIM = 100*(pathvalue_IM[:T_max]-ssvalue_IM)*scaleval / ssvalue_IM  
        dEX = 100*(pathvalue_EX[:T_max]-ssvalue_EX)*scaleval / ssvalue_EX      
        dNX = dEX-dIM  
        ax.plot(np.arange(T_max),dNX,label=label,linestyle=lstyle,color=color,linewidth=lwidth)
        ax.set_ylabel('\% diff. to s.s.')
        ax.set_title(pathlabel)

    else:

        if abs(ssvalue) > 0: 
        
            ax.plot(np.arange(T_max),((dpathvalue[:T_max])*scaleval/ssvalue)*100,label=label, linestyle=lstyle,color=color,linewidth=lwidth)
            ax.set_ylabel('$\%$ diff. to s.s.')
            ax.set_title(pathlabel)
        
        else:
        
            ax.plot(np.arange(T_max),((dpathvalue[:T_max])*scaleval)*100,label=label, linestyle=lstyle, color=color,linewidth=lwidth)
            ax.set_ylabel('$\%$ diff. to s.s.')
            ax.set_title(pathlabel)                           

def show_IRFs(models,paths,labels=None,
              T_max=17,scale=True,
              lwidth=1.3,lstyles=None,colors=None,palette=None,
              maxcol=5,figsize=None,
              lfsize=15,legend_window=0,compare_LP=None,CI=False,do_stds=False):
    """ plot IRFs for a list of models """


    # models: list of models
    # paths: list of path variables
    # labels: list of labels for each model

    # compare_LP: data when comparing to LP

    # a. models as list
    if type(models) is dict: models = [model for model in models.values()]
    model = models[0]
    par = model.par

    # b. inputs
    if type(paths) is str: paths = paths_defaults[paths]
    if labels is None: labels = [None]*len(models)
    assert len(labels) >= len(models), f'{len(labels) = } must be same as {len(models) = }'

    if T_max is None: T_max = par.T   
    T_max_LP = 17

    # c. figure
    num = len(paths)
    nrows = num//maxcol+1
    ncols = np.fmin(num,maxcol)
    if num%maxcol == 0: nrows -= 1 

    if figsize is None:
        fig = plt.figure(figsize=(4.3*ncols,3.6*nrows))
    else:
        fig = plt.figure(figsize=(figsize[0]*ncols,figsize[1]*nrows))

    for i,pathname in enumerate(paths):

        ax = fig.add_subplot(nrows,ncols,i+1)

        # axis
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        if T_max < 100:
            ax.set_xticks(np.arange(0,T_max,4))
        
        # models
        for j, model_ in enumerate(models):  
            
            if lstyles is None:
                lstyle = '-' 
            else:
                lstyle = lstyles[j]

            if colors is None:               
                color = None
            elif palette is not None:
                color = color_palette(palette,1+len(models))[1+j]
                if colors is not None:
                    color = colors[j]
            else:
                color = colors[j]

            label = labels[j]

            if j == 0: ax.plot(np.zeros(T_max), '-', color='black')
            _plot_IRFs(ax,model_,pathname,scale,lstyle,color,lwidth,label,T_max)

            ax.set_xlim([0,T_max-1])
            if i >= ncols*(nrows-1): ax.set_xlabel('Quarters',fontsize=16)
            
        if compare_LP is not None:

            try:

                ax.plot(np.arange(T_max_LP),compare_LP['IRF'][pathname][:T_max_LP]*100,label='LP', linestyle='--',color='black',linewidth=lwidth)
                if do_stds:
                    if CI:
                        prob=0.05
                        sign = norm.ppf(1-prob/2)
                        SE = compare_LP['SE'][pathname]
                        LO = compare_LP['IRF'][pathname] - SE
                        HI = compare_LP['IRF'][pathname] + SE
                        ax.fill_between(np.arange(T_max_LP), LO[:T_max_LP], HI[:T_max_LP], alpha=0.15, color='C0')
                        LO = compare_LP['IRF'][pathname] - sign*SE
                        HI = compare_LP['IRF'][pathname] + sign*SE
                        ax.fill_between(np.arange(T_max_LP), LO[:T_max_LP], HI[:T_max_LP], alpha=0.08, color='C0')
                else:
                    if CI:
                        LO = compare_LP['IRF'][pathname] + compare_LP['LO'][pathname][:,0]*100  
                        HI = compare_LP['IRF'][pathname] + compare_LP['HI'][pathname][:,0]*100  
                        ax.fill_between(np.arange(T_max_LP), LO[:T_max_LP], HI[:T_max_LP], alpha=0.25, color='C0')
                        LO = compare_LP['IRF'][pathname] + compare_LP['LO'][pathname][:,1]*100  
                        HI = compare_LP['IRF'][pathname] + compare_LP['HI'][pathname][:,1]*100  
                        ax.fill_between(np.arange(T_max_LP), LO[:T_max_LP], HI[:T_max_LP], alpha=0.15, color='C0')
            except:
                pass

        if len(models) > 1:
            if i == legend_window: ax.legend(frameon=True, prop={'size': lfsize})

    fig.tight_layout(pad=1.6)
    plt.show()

    return fig

def show_IRFs_vs_data(models,data,paths=None,labels=None,T_max=17,lwidth=2.5,lstyles=None,colors=None,maxcol=5,figsize=None,legend_window=0,filename=None):

    if paths is None: paths = paths_defaults['standard_vs_data']

    if colors is None: colors = ['royalblue', 'orange', 'darkgreen']
    if labels is None: labels = [model.name for model in models]
    if lstyles is None: lstyles = ['-','-.','--',':']
    
    fig = show_IRFs(models=models,paths=paths,labels=labels,
                    T_max=T_max,scale=True, 
                    lwidth=lwidth,lstyles=lstyles,colors=colors,maxcol=maxcol,figsize=figsize,
                    compare_LP=data.LP_IRFS_SOE,CI=True,legend_window=legend_window,do_stds=False)

    if filename is not None: fig.savefig(f'plots/{filename}.pdf')
    return fig

def show_IRFs_robust(modelspecs,paths,labels=None,
              T_max=17,scale=True,
              lwidth=1.3,lstyles=None,colors=None,palette=None,
              maxcol=4,figsize=None,
              lfsize=15,legend_window=0):

    # a. inputs
    if labels is None: labels = [None]
    if T_max is None: T_max = 17

    # c. figure
    num = len(paths)*len(modelspecs)
    nrows = num//maxcol+1
    ncols = np.fmin(num,maxcol)
    if num%maxcol == 0: nrows -= 1 

    if figsize is None:
        fig = plt.figure(figsize=(4.3*ncols,3.6*nrows))
    else:
        fig = plt.figure(figsize=(figsize[0]*ncols,figsize[1]*nrows))

    i = 0
    for sup_ylabel,models in modelspecs.items():        
        j = 0
        for pathname in paths:
            
            ax = fig.add_subplot(nrows,ncols,i+1)

            if j == 0:
            
                pad=0.8
                ax.annotate(sup_ylabel, xy=(-0.5, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size=18, ha='right', va='center', rotation=90)            
            
            for j,model_ in enumerate(models):

                # axis
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
                if T_max < 100:
                    ax.set_xticks(np.arange(0,T_max,4))
        
                if lstyles is None:
                    lstyle = '-' 
                else:
                    lstyle = lstyles[j]

                if colors is None:               
                    color = None
                elif palette is not None:
                    color = color_palette(palette,1+len(models))[1+j]
                    if colors is not None:
                        color = colors[j]
                else:
                    color = colors[j]

                label = labels[j]

                if j == 0: ax.plot(np.zeros(T_max), '-', color='black')
                _plot_IRFs(ax,model_,pathname,scale,lstyle,color,lwidth,label,T_max)

            ax.set_xlim([0,T_max-1])

            if i >= ncols*(nrows-1): ax.set_xlabel('Quarters',fontsize=16)
            if i == legend_window: ax.legend(frameon=True,prop={'size':lfsize})
            i += 1
            j += 1

    fig.tight_layout(pad=1.6)
    plt.show()

    return fig    