# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import print_function, division
import numpy as np
from scipy.interpolate import interp1d 
import warnings
import matplotlib.pyplot as plt # just for debugging
import sys # just for debugging

def ccm_dered(wave, flux, EBV, R_V=3.1, model='ccm89'):
    """Deredden a flux array using the CCM parameterization.
    
    Parameters
    ----------
    wave: np.array
        wavelength in Angstroms
    flux: np.array
    EBV: float
        E(B-V) differential extinction
    R_V: float, optional
        defaults to standard Milky Way average of 3.1
    model: string, optional
        'ccm89' is the default Cardelli, Clayton, & Mathis (1989), but does
            include the O'Donnell (1994) parameters to match IDL astrolib.
        'gcc09' is Gordon, Cartledge, & Clayton (2009).
    
    Returns
    ----------
    dered_flux: np.array
        dereddened flux vector
        
    Notes
    ----------
    Cardelli, Clayton, & Mathis (1989) parameterization is used.
    In the optical range, the updated parameters of O'Donnell (1994) are used 
        instead of CCM.
    Function valid between 910 A and 3.3 microns, although note the original 
        CCM values were derived using only >1250 A data.
    Gordon, Cartledge, & Clayton (2009) has updated UV coefficients, and is
        valid from 910 A to 3030 A. This function will use CCM89 at longer
        wavelengths if GCC09 is selected, but note that the two do not connect
        perfectly smoothly. There is a small discontinuity at 3030 A. Note that
        GCC09 equations 14 and 15 apply to all x>5.9 (not limited to x<8 as stated; K. Gordon, priv. comm.).

    References
    ----------
    Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    O'Donnell, J. E. 1994, ApJ, 422, 158O

    """
    
    model = model.lower()
    if model not in ['ccm89','gcc09']:
        raise ValueError('ccm_dered: model must be ccm89 or gcc09')
    
    x = 1e4 / wave      # inverse microns
    
    if any(x < 0.3) or any(x > 11):
        raise ValueError('ccm_dered valid only for wavelengths from 910 A to '+
            '3.3 microns')
    if any(x > 8) and (model == 'ccm89'):
        warnings.warn('CCM89 should not be used below 1250 A.')
#    if any(x < 3.3) and any(x > 3.3) and (model == 'gcc09'):
#        warnings.warn('GCC09 has a discontinuity at 3030 A.')
    
    a = np.zeros(x.size)
    b = np.zeros(x.size)
    
    # NIR
    valid = np.where((0.3 <= x) & (x < 1.1))
    a[valid] =  0.574 * x[valid]**1.61
    b[valid] = -0.527 * x[valid]**1.61
    
    # optical, using O'Donnell (1994) values
    valid = np.where((1.1 <= x) & (x < 3.3))
    y = x[valid] - 1.82
    coef_a = np.array([-0.505, 1.647, -0.827, -1.718, 1.137, 0.701, -0.609,
        0.104, 1.])
    coef_b = np.array([3.347, -10.805, 5.491, 11.102, -7.985, -3.989, 2.908,
        1.952, 0.])
    function_a = np.poly1d(coef_a)
    function_b = np.poly1d(coef_b)
    a[valid] = function_a(y)
    b[valid] = function_b(y)
    
    # UV
    valid = np.where((3.3 <= x) & (x < 8))
    y = x[valid]
    F_a = np.zeros(valid[0].size)
    F_b = np.zeros(valid[0].size)
    select = np.where(y >= 5.9)
    yselect = y[select] - 5.9
    F_a[select] = -0.04473 * yselect**2 - 0.009779 * yselect**3
    F_b[select] = 0.2130 * yselect**2 + 0.1207 * yselect**3
    a[valid] = 1.752 - 0.316*y - (0.104 / ((y-4.67)**2 + 0.341)) + F_a
    b[valid] = -3.090 + 1.825*y + (1.206 / ((y-4.62)**2 + 0.263)) + F_b
    
    # far-UV CCM89 extrapolation
    valid = np.where((8 <= x) & (x < 11))
    y = x[valid] - 8.
    coef_a = np.array([-0.070, 0.137, -0.628, -1.073])
    coef_b = np.array([0.374, -0.420, 4.257, 13.670])  
    function_a = np.poly1d(coef_a)
    function_b = np.poly1d(coef_b)
    a[valid] = function_a(y)
    b[valid] = function_b(y)
    
    # Overwrite UV with GCC09 model if applicable. Not an extrapolation.
    if model == 'gcc09':
        valid = np.where((3.3 <= x) & (x < 11))
        y = x[valid]
        F_a = np.zeros(valid[0].size)
        F_b = np.zeros(valid[0].size)
        select = np.where(5.9 <= y)
        yselect = y[select] - 5.9
        F_a[select] = -0.110 * yselect**2 - 0.0099 * yselect**3
        F_b[select] = 0.537 * yselect**2 + 0.0530 * yselect**3
        a[valid] = 1.896 - 0.372*y - (0.0108 / ((y-4.57)**2 + 0.0422)) + F_a
        b[valid] = -3.503 + 2.057*y + (0.718 / ((y-4.59)**2 + 0.0530)) + F_b
    
    
    A_V = EBV * R_V
    A_lambda = A_V * (a + b/R_V)
    dered_flux = flux * 10**(0.4 * A_lambda)
    
#    return dered_flux
    return A_lambda / A_V
    
    
def fm_dered(wave, flux, EBV, R_V=3.1, model='f99'):
    """Deredden a flux array using the Fitzpatrick & Massa parameterization.
    
    Parameters
    ----------
    wave: np.array
        wavelength in Angstroms
    flux: np.array
    EBV: float
        E(B-V) extinction
    R_V: float, optional
        defaults to standard Milky Way average of 3.1
    model: string, optional
        'f99' is the default Fitzpatrick (1999)
        'fm07' is Fitzpatrick & Massa (2007)
    
    Returns
    ----------
    dered_flux: np.array
            dereddened flux vector
        
    Notes
    ----------
    Uses Fitzpatrick (1999) by default, which relies on the UV parametrization of
    Fitzpatrick & Massa (1990) and spline fitting in the optical and IR. This function is defined from 
    910 A to 6 microns, but note the claimed validity goes down only to 
    1150 A.
        
    The fm07 model uses the Fitzpatrick & Massa (2007) parametrization,
    which has a slightly different functional form. That paper claims it
    preferable, although it is unclear if signficantly (Gordon et al. 2009). It is not the
    literature standard, so not default here.
    
    The gcc09 model uses the same functional form as f99, but updates the UV
    parameters c4 and gamma using FUSE data for the <1150 A region.

    References
    ----------
    Fitzpatrick, E. L. 1999, PASP, 111, 63
    Fitpatrick, E. L. & Massa, D. 1990, ApJS, 72, 163
    Fitpatrick, E. L. & Massa, D. 2007
    Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320

    """
    
    model = model.lower()
    if model not in ['f99','fm07']:
        raise ValueError('fm_dered: model must be f99 or fm07')
    
    x = 1e4 / wave      # inverse microns
    
    if any(x < 0.167) or any(x > 11):
        raise ValueError('fm_dered valid only for wavelengths from 910 A to '+
            '6 microns')
    
    # UV region
    uvsplit = 10000. / 2700.  # Turn 2700A split into inverse microns.
    uv_region = np.where(x >= uvsplit)
    y = x[uv_region]
    k_uv = np.zeros(y.size)
    
    # Fitzpatrick (1999) model
    if model in ['f99','gcc09']:
        x0, gamma = 4.596, 0.99
        c3, c4 = 3.23, 0.41
        if model == 'gcc09':
            c4 = 0.377 # update with real value
            gamma = 0.99 # update with real value
        c2 = -0.824 + 4.717 / R_V
        c1 = 2.030 - 3.007 * c2
        D = y**2 / ((y**2-x0**2)**2 + y**2 * gamma**2)
        F = np.zeros(y.size)
        valid = np.where(y >= 5.9)
        F[valid] = 0.5392 * (y[valid]-5.9)**2 + 0.05644 * (y[valid]-5.9)**3
        k_uv = c1 + c2*y + c3*D + c4*F
    # Fitzpatrick & Massa (2007) model
    if model == 'fm07':
        x0, gamma = 4.592, 0.922
        c1, c2, c3, c4, c5 = -0.175, 0.807, 2.991, 0.319, 6.097
        D = y**2 / ((y**2-x0**2)**2 + y**2 * gamma**2)
        valid = np.where(y <= c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid]
        valid = np.where(y > c5)
        k_uv[valid] = c1 + c2*y[valid] + c3*D[valid] + c4*(y[valid]-c5)**2
        
#        plt.plot(y,k_uv) # debug
#        plt.show() # debug
#        sys.exit() # debug
        
    # Calculate values for UV spline points to anchor OIR fit
    x_uv_spline = 10000. / np.array([2700., 2600.])
    D = x_uv_spline**2 / ((x_uv_spline**2-x0**2)**2 + x_uv_spline**2 * gamma**2)
    k_uv_spline = c1 + c2*x_uv_spline +c3*D
    print('k_uv_spline',k_uv_spline) #debug
        
    # Optical / IR
    OIR_region = np.where(x < uvsplit)
    y = x[OIR_region]
    k_OIR = np.zeros(y.size)    
    
    # Fitzpatrick (1999) model
    if model == 'f99':
#        anchors_extinction = np.array([0, 0.265*R_V/3.1, 0.829*R_V/3.1,  # IR
#            -0.426 + 1.0044*R_V,  # optical
#            -0.050 + 1.0016*R_V,
#            0.701 + 1.0016*R_V,
#            -1.208 + 1.0032*R_V - 0.00033*R_V**2]) # this equation is wrong in F99 - it does not yield the correct value
        anchors_extinction = np.array([0, 0.26469*R_V/3.1, 0.82925*R_V/3.1,  # IR # these are the IDL astrolib values, not the F99 values
            -0.422809 + 1.00270*R_V + 2.13572e-04**R_V**2, # optical
            -5.13540e-02 + 1.00216*R_V - 7.35778e-05*R_V**2,
            0.700127 + 1.00184*R_V - 3.32598e-05*R_V**2,
            (1.19456 + 1.01707*R_V - 5.46959e-03*R_V**2 + 7.97809e-04*R_V**3 + 
                -4.45636e-05**R_V**4)])
        anchors_k = np.append(anchors_extinction-R_V, k_uv_spline)
        print('anchors_k+R_V',anchors_k+R_V)  # debug
        anchors_x = 1e4 / np.array([26500., 12200., 6000., 5470., 4670., 4110.])
        anchors_x = np.append(0., anchors_x)  # For well-behaved spline.
        anchors_x = np.append(anchors_x, x_uv_spline)
        OIR_spline = interp1d(anchors_x, anchors_k, kind='cubic') 
        k_OIR = OIR_spline(y) 
    # Fitzpatrick & Massa (2007) model
    if model == 'fm07':
        anchors_k_opt = np.array([0., 1.322, 2.055])
        IR_wave = np.array([float('inf'), 4., 2., 1.333, 1.])
        anchors_k_IR = (-0.83 + 0.63*R_V) * IR_wave**-1.84 - R_V
        anchors_k = np.append(anchors_k_IR, anchors_k_opt)
        anchors_k = np.append(anchors_k, k_uv_spline)
        anchors_x = np.array([0., 0.25, 0.50, 0.75, 1.])  # IR
        opt_x = 1e4 / np.array([5530., 4000., 3300.])  # optical
        anchors_x = np.append(anchors_x, opt_x)
        anchors_x = np.append(anchors_x, x_uv_spline)
        
        print('anchors_x=',anchors_x) # debug
        print('anchors_k=',anchors_k) # debug
        
        OIR_spline = interp1d(anchors_x, anchors_k, kind='cubic') 
        k_OIR = OIR_spline(y) 
        
    k = np.append(k_uv, k_OIR)
    
    dered_flux = flux * 10**(0.4 * EBV * (k+R_V))
    
#    return dered_flux
    return (k+R_V) / R_V    # A_lambda / A_V