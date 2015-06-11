
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from IPython.html.widgets import interact, fixed

from solve_system import derivs, solve_system, get_sol
from initial_conditions import star, star_cond, S_cond, S_condd

def plot_pos(rx, ry, sx, sy, j=0):
    """Plots the solved positions of the system at any time
    
    Parameters
    ----------
    rx: array, x postitions of all stars
    ry: array, y postitions of all stars
    sx: array, x postitions of S
    sy: array, y postitions of S
    j: float, time at which you are plotting
    
    Returns
    -------
    Scatter plot of the system at any given time.
    """
    plt.figure(figsize=(6,6))
    
    plt.scatter(0,0,color= 'c', label= 'M')
    plt.scatter(sx[j], sy[j], color = 'b', label = 'S')
    plt.scatter(rx[j], ry[j], color = 'g', marker = "*", label = 'stars', s=6)
        
    plt.legend(loc='upper left')
    plt.xlim(-65,65)
    plt.ylim(-65,65)
    plt.box(False)
    plt.tick_params(axis = 'x', top = 'off', bottom = "off", labelbottom= 'off')
    plt.tick_params(axis = 'y', right = 'off', left= "off", labelleft= 'off')
    plt.show()
    
def plot_pos_cm(rx, ry, sx, sy, M, S, j=0):
    """Plots the solved positions of the system at any time around the system's center of mass
    
    Parameters
    ----------
    rx: array, x postitions of all stars
    ry: array, y postitions of all stars
    sx: array, x postitions of S
    sy: array, y postitions of S
    M,S: float, respective mass of galactic neucli
    j: float, time at which you are plotting
    
    Returns
    -------
    Scatter plot of the system at any given time.
    """
    plt.figure(figsize=(6,6))
    
    plt.scatter(0,0,color= 'k', label= 'CM', marker = "+", s= 100)
    plt.scatter((sx[j]-((M*sx[j])/(M+S))), (sy[j]-((M*sy[j])/(M+S))), color = 'b',label ='S')
    plt.scatter((rx[j]-((S*sx[j])/(M+S))),(ry[j]-((S*sy[j])/(M+S))),color='g',marker="*", label='stars', s=6)
    plt.scatter((-(M*sx[j])/(M+S)),(-(M*sy[j])/(M+S)), color = 'c', label = "M")
        
    plt.legend(loc='upper left')
    plt.xlim(-65,65)
    plt.ylim(-65,65)
    plt.box(False)
    plt.tick_params(axis = 'x', top = 'off', bottom = "off", labelbottom= 'off')
    plt.tick_params(axis = 'y', right = 'off', left= "off", labelleft= 'off')
    plt.show()   
    
def sub_retro(rx, ry, sx, sy, times):
    """Creates a 2 by 3 subplot of the system at 6 chosen times
    
    Parameters
    ----------
    rx: array, x postitions of all stars
    ry: array, y postitions of all stars
    sx: array, x postitions of S
    sy: array, y postitions of S
    times: float, time at which you are plotting
    
    Returns
    -------
    2 by 3 subplot of the system at 6 chosen times
    """
    fig, ax  = plt.subplots(2,3, figsize=(12,8))
    c=1
    for i in range(2):
        for j in range(3):
            plt.sca(ax[i,j])
            t=times[i,j]
            plt.scatter(0,0,color= 'c', label= 'M', s=40)
            plt.scatter(sx[t], sy[t], color = 'b', label = 'S', s=40)
            plt.scatter(rx[t], ry[t], color = 'g', label = 'stars', s=7)
            plt.xlabel(c+2, fontsize=15)
            plt.xlim(-45,45)
            plt.ylim(-45,45)
            plt.box(False)
            plt.tick_params(axis = 'x', top = 'off', bottom = "off", labelbottom= 'off')
            plt.tick_params(axis = 'y', right = 'off', left= "off", labelleft= 'off')
            c+=1
    plt.tight_layout()
    
    
def sub_direct(rx, ry, sx, sy, time):
    """Creates a 2 by 4 subplot of the system at 7 chosen times
    
    Parameters
    ----------
    rx: array, x postitions of all stars
    ry: array, y postitions of all stars
    sx: array, x postitions of S
    sy: array, y postitions of S
    times: float, time at which you are plotting
    
    Returns
    -------
    2 by 4 subplot of the system at 7 chosen times
    """
    c=1
    fig, ax = plt.subplots(2,4, figsize=(20,10))
    for i in range(2):
        for j in range(4):
            plt.sca(ax[i,j])
            t = time[i,j]
            plt.scatter(0,0,color= 'b', label= 'M', s=30)
            plt.scatter(sx[t], sy[t], color = 'r', label = 'S', s=30)
            plt.scatter(rx[t], ry[t], color = 'g', label = 'stars', s=5)
            plt.xlabel(c-3, fontsize=15)
            plt.xlim(-60,60)
            plt.ylim(-60,60)
            plt.box(False)
            plt.tick_params(axis = 'x', top = 'off', bottom = "off", labelbottom= 'off')
            plt.tick_params(axis = 'y', right = 'off', left= "off", labelleft= 'off')
            c+=1
    plt.tight_layout()