def buckley_leverett1d(nt, Sw0, L, nx, sigma, bc_value, muw, muo, q, A, poro, Sw_data, krw_data, kro_data):
    """ 
    Solve Buckley-Leverett PDE using forward-time/backward-space scheme
    """
    import numpy
    import matplotlib.pyplot as pyplot

    def interstitial_velocity(q, A, poro):
        # interstitial velocity vt
        vt = q / A * poro
        return vt

    def fractional_flow(krw, muw, kro, muo):
        # fractional flow Fww
        Fww = 1 / (1 + (kro / muo) * (muw / krw))
        return Fww    
    
    # calculate interstitial velocity
    vt = interstitial_velocity(q, A, poro)
    
    # calculate dx
    dx = L / (nx - 1)
    
    # calculate dt from CFwL
    dt = sigma * dx / vt  # time-step size 
    
    # Discretize the domain.
    x = numpy.linspace(0.0, L, num=nx)   
    
    # integrate solution in time
    Sw_hist = [Sw0.copy()]
    Sw = Sw0.copy()
    for n in range(nt):
        # Compute the fractional flow.
        krw, kro = interpolate_relperm(Sw_data, krw_data, kro_data, Sw)
        Fw = fractional_flow(krw, muw, kro, muo)
        
        # Advance in time.
        Sw[1:] = Sw[1:] - (vt * dt / dx) * (Fw[1:] - Fw[:-1])
        
        # Set the left boundary condition.
        Sw[0] = bc_value
        
        # Record the time-step solution.
        Sw_hist.append(Sw.copy()) 

    # Plot Sw over x
    # fig = pyplot.figure(figsize=(6.0, 4.0))
    pyplot.xlabel(r'$x$')
    pyplot.ylabel(r'Sw')
    pyplot.grid()
    pyplot.plot(x, Sw0, label='Initial',
                color='C0', linestyle='--', linewidth=2)
    pyplot.plot(x, Sw, label='nt = {}'.format(nt),
                color='C1', linestyle='-', linewidth=2)
    pyplot.xlim(0.0, L)
    pyplot.ylim(0, 1)  
    pyplot.legend()
    pyplot.grid()
