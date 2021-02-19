def interpolate_relperm(Sw_data, krw_data, kro_data, Sw_new):
    """
    Spline interpolation of relative permeability data
    """
    from scipy import interpolate

    # Spline interpolation of data
    krw_interp = interpolate.splrep(Sw_data, krw_data, s=0)
    kro_interp = interpolate.splrep(Sw_data, kro_data, s=0)    

    # Interpolate krw and kro at given Sw
    krw_new = interpolate.splev(Sw_new, krw_interp, der=0)
    kro_new = interpolate.splev(Sw_new, kro_interp, der=0)

    return krw_new, kro_new
