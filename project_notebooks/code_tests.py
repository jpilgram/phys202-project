def E_cons(rS, t, N, M, S):
    """Tests energy of system at any point
    
    Parameters
    ----------
    rS: float, initial radius of S
    t: array, array of times
    N: int, number of stars in each layer
    M: float, mass of M
    S: float, mass of S
    
    Returns
    -------
    Returns total energy for each timestep"""
    
    G = 4.498e-6
    m=1
    
    sx,sy,svx,svy = S_cond(rS,45)
    c = star_cond(N, M, 0.6*rS)
    
    for i in range(N):
        ic = c[i]
        ic = np.append(ic, np.array([sx, svx, sy ,svy]))

        sol = solve_system(ic, t, M, S)
        rx = np.array([i[0] for i in sol])
        ry = np.array([i[2] for i in sol])
        rdx = np.array([i[1] for i in sol])
        rdy = np.array([i[3] for i in sol])
        Rx = np.array([i[4] for i in sol])
        Ry = np.array([i[6] for i in sol])
        Rdx = np.array([i[5] for i in sol])
        Rdy = np.array([i[7] for i in sol])
        rhox = Rx - rx
        rhoy = Ry -ry
        U= -(G*M*S)/(np.sqrt(Rx**2+Ry**2)+ (-G*M)/(np.sqrt(rx**2+ry**2)+ -(G*S)/np.sqrt(rhox**2+rhoy**2))) 
        vS = np.sqrt((Rdx)**2+(Rdy**2))
        vm = np.sqrt((rdx)**2+(rdy)**2)
        K  = 0.5*S*vS**2 + 0.5*m*vm**2
    return K+U

