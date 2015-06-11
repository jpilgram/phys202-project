import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed

from initial_conditions import star, star_cond, S_cond, S_condd

def derivs(rvec, t, M, S):
    """Computes the derivs for the 3 body system at rvec(t)
    Parameters
    ----------
    rvec: np.array, solution vector
    M, S: float, central mass of respective galaxies
    t: array, times at which want derivatives
    
    Returns
    -------
    vector: the derivatives of function
    """
    
    rx = rvec[0]
    ry = rvec[2]
    rdx = rvec[1]
    rdy = rvec[3]
    
    Rx = rvec[4]
    Ry = rvec[6]
    Rdx = rvec[5]
    Rdy = rvec[7]
    
    rhox = Rx - rx
    rhoy = Ry - ry
    
    R_mag= np.sqrt(Rx**2 + Ry**2)
    r_mag = np.sqrt(rx**2 + ry**2)
    rho_mag = np.sqrt(rhox**2 + rhoy**2)
    
    G = 4.498e-6
    
    rdx2 = -G*((M/r_mag**3)*rx - rhox*(S/rho_mag**3)  + Rx*(S/R_mag**3))
    rdy2 = -G*((M/r_mag**3)*ry - rhoy*(S/rho_mag**3)  + Ry*(S/R_mag**3))
    
    Rdx2 = -G*(Rx*((M+S)/R_mag**3))
    Rdy2 = -G*(Ry*((M+S)/R_mag**3))
    
    return np.array([rdx, rdx2, rdy, rdy2, Rdx, Rdx2, Rdy, Rdy2])

def solve_system(ic, t, M, S):
    """Computes solution of 3 Body Problem for single initial condition ic
    Parameters
    ----------
    ic: array or list, initial conditions of system
    t: array, array of times for which the system will be solved
    M, S: float, central mass of respective galaxies
    Returns
    -------
    returns array of the solution vectors of the differential equation
    """
    
    return odeint(derivs, ic, t, args=(M, S), atol=1e-4, rtol=1e-5)

def get_sol(rS, t, M, S, N):
    """Solves the system for N number of stars
    
    Parameters
    ----------
    rS: float, initial radius of S
    t: array, array of times
    M: float, mass of M
    S: float, mass of S
    N: int, number of stars wanted in inner layer,
        each consecutive layer will have 6 additional stars
    
    Returns
    -------
    Returns the solutions postitions for each star and S at all times for a retrograde passage
    """
    G = 4.498e-6

    sx,svx,sy,svy = S_cond(rS,50, M, S)

    radi = np.array([0.2*rS,0.3*rS,0.4*rS,0.5*rS,0.6*rS])
    k=N+1
    star_sols = []
    for r in radi:
        c = star_cond(k,M,r)
        for i in range(k):
            ic = c[i]
            ic = np.append(ic, np.array([sx, svx, sy, svy]))
            sol = solve_system(ic, t, M, S)

            x = sol[:,0]
            y = sol[:,2]
            X = sol[:,4]
            Y = sol[:,6]
            star_sols.append(x)
            star_sols.append(y)
        k+=6
    all_stars = np.transpose(np.vstack(star_sols))
    starx = all_stars[:,0::2]
    stary = all_stars[:,1::2]
        
    return starx, stary, X, Y

def get_sol_direct(rS, t, M, S, N):
    """Solves for a direct passage of the system for N nummber of stars for inner layer
    Parameters
    ----------
    rS: float, initial radius of S
    t: array, array of times
    M: float, mass of M
    S: float, mass of S
    N: int, number of stars in each layer inner layer, 6 more stars in each additional layer
    
    Returns
    -------
    returns the solution postitions for all the stars and S at all times for a direct passage
    """
    
    G = 4.498e-6

    sx,svx,sy,svy = S_condd(rS,50, M, S)

    radi = np.array([0.2*rS,0.3*rS,0.4*rS,0.5*rS,0.6*rS])
    k=N+1
    star_sols = []
    for r in radi:
        c = star_cond(k,M,r)
        for i in range(k):
            ic = c[i]
            ic = np.append(ic, np.array([sx, svx, sy, svy]))
            sol = solve_system(ic, t, M, S)

            x = sol[:,0]
            y = sol[:,2]
            X = sol[:,4]
            Y = sol[:,6]
            star_sols.append(x)
            star_sols.append(y)
        k+=6
    all_stars = np.transpose(np.vstack(star_sols))
    starx = all_stars[:,0::2]
    stary = all_stars[:,1::2]
        
    return starx, stary, X, Y
