import time
import numpy as np

import matplotlib.pyplot as plt
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.rcParams.update({"axes.grid" : True, "grid.color": "black", "grid.alpha":"0.25", "grid.linestyle": "--"})
plt.rcParams.update({'font.size': 12})

from consav import elapsed

def create_fig(figsize=(6,6/1.5)):
    
    fig = plt.figure(figsize=figsize,dpi=100)
    ax = fig.add_subplot(1,1,1)

    return fig,ax

def save_fig(fig,ax,filename=None,title=None,legend=False,ncol=1,bbox=False,ylabel='',T_max=48):

    if legend: ax.legend(frameon=True)
    ax.set_xlabel('months')
    ax.set_ylabel(ylabel)
    ax.set_xticks(np.arange(T_max+1)[::12])
    ax.set_xlim([0,T_max]);
    
    if not filename is None:
        
        fig.tight_layout()
        fig.savefig(f'results/{filename}.pdf')
    
        if not title is None:
            
            ax.set_title(title,fontsize=14,weight='bold')

            if legend: 
                if bbox:                
                    ax.legend(frameon=True,ncol=ncol,fontsize=12,bbox_to_anchor=(1.05,0.4))
                else:
                    ax.legend(frameon=True,ncol=ncol,fontsize=12)
                    
            fig.tight_layout()
            fig.savefig(f'results/{filename}_titled.pdf')    

def vary_par(model,parname,values,do_ss=True,skip_hh=False,fac=False,shocks=None,
    do_print=False,show_info=False,do_print_full=False):
    
    t0 = time.time()
    if do_print: print(f'baseline: {parname} = {model.par.__dict__[parname]:.4f}')
    
    models = []
    for value in values:
    
        if do_print: print(f'{parname} = {value:7.4f}',end='')

        try:

            model_ = model.copy()
            
            if fac:
                model_.par.__dict__[parname] *= value
            else:
                model_.par.__dict__[parname] = value
            
            if do_ss: model_.find_ss(calib_beta=False,do_print=do_print_full)

            model_.compute_jacs(do_print=do_print_full,skip_hh=skip_hh,skip_shocks=True)
            model_.find_transition_path(do_print=do_print_full,shocks=shocks,do_end_check=False)    
            model_.calc_moms_path(do_print=do_print_full)
            
            if do_print and show_info:
                print(f': var_u = {model_.moms["var_u"]:5.2f}')
            else:
                print('')

            models.append(model_)

        except Exception as e:

            print(f': error, {e}')
            
    if do_print: print(f'completed in {elapsed(t0)}')
            
    return models
  
def IRF_figs(models,labels,lss,colors,title,prefix,parname,varname,legend=False,T_max=48,ylim=None):
    
    if lss is None: lss = ['-']*len(models)
    if colors is None: colors = ['black']*len(models)
     
    fig,ax = create_fig()

    for model,label,ls,color in zip(models,labels,lss,colors):

        IRF,ylabel = model.get_IRF(varname)
        ax.plot(IRF,ls=ls,color=color,label=label)
    
    if not ylim is None: ax.set_ylim(ylim)
    save_fig(fig,ax,f'{prefix}_{varname}_{parname}',
        title=title,ylabel=ylabel,legend=legend,T_max=T_max)
   