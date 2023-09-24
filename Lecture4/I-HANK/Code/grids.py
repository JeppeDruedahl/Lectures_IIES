import numpy as np

##########
# income # 
##########

def variance(x, pi):
    """ Variance of discretized random variable with support x and probability mass function pi."""

    return np.sum(pi * (x - np.sum(pi * x)) ** 2)
    
def markov_rouwenhorst(rho, sigma, N=7):
    """ Rouwenhorst method analog to markov_tauchen"""

    # parametrize Rouwenhorst for n=2
    p = (1 + rho) / 2
    Pi = np.array([[p, 1 - p], [1 - p, p]])

    # implement recursion to build from n=3 to n=N
    for n in range(3, N + 1):
        P1, P2, P3, P4 = (np.zeros((n, n)) for _ in range(4))
        P1[:-1, :-1] = p * Pi
        P2[:-1, 1:] = (1 - p) * Pi
        P3[1:, :-1] = (1 - p) * Pi
        P4[1:, 1:] = p * Pi
        Pi = P1 + P2 + P3 + P4
        Pi[1:-1] /= 2

    # invariant distribution and scaling
    pi = stationary(Pi)
    s = np.linspace(-1, 1, N)
    s *= (sigma     / np.sqrt(variance(s, pi)))
    y = np.exp(s) / np.sum(pi * np.exp(s))

    return y, pi, Pi

def stationary(Pi, pi_seed=None, tol=1E-11, maxit=10_000):
    """Find invariant distribution of a Markov chain by iteration."""

    if pi_seed is None:
        pi = np.ones(Pi.shape[0]) / Pi.shape[0]
    else:
        pi = pi_seed

    for it in range(maxit):
        pi_new = pi @ Pi
        if np.max(np.abs(pi_new - pi)) < tol:
            break
        pi = pi_new
    else:
        raise ValueError(f'No convergence after {maxit} forward iterations!')
    
    pi = pi_new

    return pi
    
#########
# grids #
#########

def asset_grid(amin, amax, n):
    """ create asset grid """

    # find maximum ubar of uniform grid corresponding to desired maximum amax of asset grid
    ubar = np.log(1 + np.log(1 + amax - amin))
    
    # make uniform grid
    u_grid = np.linspace(0, ubar, n)
    
    # double-exponentiate uniform grid and add amin to get grid from amin to amax
    return amin + np.exp(np.exp(u_grid) - 1) - 1

def create_grids_full(model):
    """ create grids """

    # note: only fills out already allocated arrays

    par = model.par
    ss = model.ss
    
    # a. beta
    par.beta_grid[:] = np.linspace(par.beta_mean-par.beta_delta,par.beta_mean+par.beta_delta,par.Nbeta)
        
    # b. a
    par.a_grid[:] = asset_grid(par.a_min,par.a_max,par.Na)
    assert par.a_min >= 0, 'a_min must be non-negative'

    if model.par.HH_type == 'HA-2A': par.ai_grid[:] = np.linspace(par.ai_min,par.ai_max,par.Nai)

    # c. earnings 
    par.e_grid_ss[:],par.e_ergodic_ss[:],par.e_trans_ss[:,:] = markov_rouwenhorst(par.rho_e,par.sigma_e,par.Ne)
     
    # d. sectors 
    par.s_set[:par.Ne] = 0 
    par.s_set[par.Ne:] = 1 
    
    par.sT_vec = np.array([0.5+par.sTransferT,1.5-par.sTransferT]) # sum to one
    
    par.s_ergodic_ss[0] = par.sT
    par.s_ergodic_ss[1] = 1-par.sT
    par.s_trans_ss = np.identity(par.Ns)

    # e. construct z_trans
    for i_s in range(par.Ns):
        par.z_ergodic_ss[i_s*par.Ne:(1+i_s)*par.Ne] = par.e_ergodic_ss * par.s_ergodic_ss[i_s]
        ss.z_trans[:,:,:] = np.kron(par.s_trans_ss, par.e_trans_ss)[np.newaxis,:,:]  
    
    # f. construct z_grid and tau_grid
    for i_s in range(par.Ns):
        if i_s == 0:
            par.z_grid_ss[i_s*par.Ne:(1+i_s)*par.Ne] = par.e_grid_ss   
        elif i_s == 1:
            par.z_grid_ss[i_s*par.Ne:(1+i_s)*par.Ne] = par.e_grid_ss  

    for i_z in range(par.Nz):
        par.tau_grid[i_z] = 1 / np.sum(par.z_ergodic_ss*par.z_grid_ss)*par.z_grid_ss[i_z]