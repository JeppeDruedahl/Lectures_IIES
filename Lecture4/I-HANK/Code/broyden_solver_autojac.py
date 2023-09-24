import numpy as np

def broyden_solver_autojac(f,x0,tol=1E-10,maxiter=500,backtrack_c=0.5,h_jac=1E-6):
    """Solve f(x)=0 using Broyden's method with automatic Jacobian calculation."""

    x = x0
    y = f(x)

    # a. initialize
    nx = x.shape[0]
    ny = y.shape[0]
    J = np.empty((nx, ny))

    for i in range(nx):
        dx = h_jac * (np.arange(nx) == i)
        J[:, i] = (f(x + dx) - y) / h_jac

    # b. iterate
    for it in range(maxiter):
        
        if np.max(np.abs(y)) < tol: return x, y
        
        # i. newton step
        dx = np.linalg.solve(J, -y)

        # ii. update
        for _ in range(30):

            try:
                ynew = f(x + dx)
            except ValueError:
                dx *= backtrack_c
            else:
                dy = ynew-y
                J = J + np.outer(((dy - J @ dx) / np.linalg.norm(dx) ** 2), dx)
                y = ynew
                x += dx
                break
            
        else:

            raise ValueError('Too many backtracks, maybe bad initial guess?')

    else:

        raise ValueError(f'No convergence after {maxiter} iterations')