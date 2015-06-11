import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed


def star(r, M, theta):
    """Computes the inital conditions of a single star
    
    Parameters
    ----------
    r: float, radius of star
    M: float, Mass of M
    theta: float, angle from +x axis at which star intially starts
    
    Returns
    -------
    returns a list of the initial conditions for a single star
    """
    
    G = 4.498e-6
    
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    
    v = np.sqrt((G*M)/np.sqrt(x**2+y**2))
    
    vx = -v*np.sin(theta)
    vy = v*np.cos(theta)
    
    return [x,vx,y,vy]

def star_cond(N, M, r):
    """Creates initial conditions for number of stars you want in inner layer
   
   Parameters
    ----------
    N: int, number of stars wanted in each shell
    M: float, mass of M
    r: radius at which you want the stars
    
    Returns
    -------
    returns array consisting of the inital conditions for each star
    """
    theta = np.linspace(0,2*np.pi,N)
    cond = np.array([star(r, M, theta[i]) for i in range(N)])
    return cond

def S_cond(R,y, M, S):
    """Computes the initial conditions for galactic nucli S
    
    Parameters
    ----------
    R: float, minimum distance between S and M
    y: flat, starting y position for S
    M,S: float, respective mass of galactic centers
    
    Returns
    -------
    x starting postion, x velocity, y startig position, y velocity
    """
    G = 4.498e-6
    x = -y**2/(4*R) + R
    theta = np.arctan(50/y)
    
    v = np.sqrt(2*(G*(M+S))/np.sqrt(x**2+y**2))
    vx = v*np.cos(theta)
    vy = -v*np.sin(theta)
    
    return x, vx, y, vy

def S_condd(R,y, M, S):
    """Computes the initial conditions for galactic nucli S
    
    Parameters
    ----------
    R: float, minimum distance between S and M
    y: flat, starting y position for S
    M, S: float: respective mass of galactic centers
    
    Returns
    -------
    x starting postion, x velocity, y startig position, y velocity
    """
    G = 4.498e-6
    x = y**2/(4*R) - R
    theta = np.arctan(50/y)
    
    v = np.sqrt(2*(G*(M+S))/np.sqrt(x**2+y**2))
    vx = -v*np.cos(theta)
    vy = -v*np.sin(theta)
    
    return x, vx, y, vy