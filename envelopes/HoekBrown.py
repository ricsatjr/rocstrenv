# =============================================================================
# DESCRIPTION:
# Module computing shear strength envelope from a given set of 
# Hoek-Brown strength parameters and a range of normal stresses. Also returns
# stress-dependent Mohr-Coulomb friction angle and cohesion values.   
# 
# REFERENCE
# Hoek, E., Carranza-Torres, C., & Corkum, B. (2002). Hoek-Brown failure 
# criterion-2002 edition. Proceedings of NARMS-Tac, 1, 267-273.
# =============================================================================

# libraries and modules
import numpy as np
from scipy.interpolate import UnivariateSpline



def HBenv(GSI=60,D=0.5,mi=15,sig_ci=25,uw=0.024,H=20):
# =============================================================================
#     Input:    
#     GSI       int or float; Geologic Strength Index (unitless)
#     D         int or float; Disturbance factor (unitless)
#     mi        int or float; material constant (unitless)
#     sig_ci    int or float; intact rock uniaxial compressive strength (MPa)
#     uw        int or float; intact rock unit weight (MN/m^3)
#     H         int or float; height of slope (m); 
#     
#     Output:    
#     sign_tau    UnivariateSpline defined by normal and shear stresses
#     sig_t       float, minimum normal stress used to define sign_tau
#     5*sig3max   float, maximum normal stress used to define sign_tau
#       
#     Description:        
#     Function returns a spline representing the Hoek-Brown failure envelope 
#     in normal stress (sign) - shear stress (tau) space, as well as the 
#     minimum and maximum normal stress used to define the spline. Normal 
#     stress range is defined from the slope height, H. 
#     
#     The spline function allows the computation of shear stresses from normal 
#     stress values not used in defining the original spline, but which are 
#     still within the range of the orgininal spline. 
#      
#     Variables and equations are from Hoek et al. (2002).
# =============================================================================
 
#    converting to float
    GSI,D,mi,sig_ci,uw,H=float(GSI),float(D),float(mi),float(sig_ci),float(uw),float(H)
    
#    1. computing Hoek-Brown derived values
    try:        
#        eq. 3
        mb=mi*np.exp((GSI-100)/(28-(14*D)))
#        eq. 4
        s=np.exp((GSI-100)/(9-(3*D)))
#        eq. 5
        a=0.5+(np.exp(-GSI/15.)-np.exp(-20/3.))/6.
#        eq. 7
        sig_t=-s*sig_ci/mb
#        eq. 17
        sig_cm=sig_ci*(mb+4*s-a*(mb-8*s))*np.power((0.25*mb+s),(a-1))/(2*(1+a)*(2+a))
#        eq. 19 (for slopes)
        sig3max=0.72*sig_cm*np.power(sig_cm/(uw*H),-0.91)
    except: 
        print "Error computing intermediate Hoek-Brown values"
#    2. creating Hoek-Brown failure envelop as sigma3 vs sigma1
    #sig3=np.sort(np.concatenate((np.linspace(sig_t*.9999,-0.9999*sig_t,1000),np.linspace(-sig_t,5*sig3max,200))))
    sig3=np.sort(np.concatenate((np.linspace(sig_t*.9999,-0.9999*sig_t,200),np.logspace(np.log(-sig_t),np.log(5*sig3max),50))))
#    eq. 2
    sig1=sig3+sig_ci*np.power((mb*sig3/sig_ci + s),a)
#    eq. 10
    ds1s3=(1+a*mb*np.power(((mb*sig3/sig_ci)+s),(a-1)))
#    3. deriving corresponding sign and tau values from sigma3 vs sigma1
#    eq. 8
    sign=((sig1+sig3)/2. - ((sig1-sig3)/2.) *(ds1s3-1)/(ds1s3+1))
#    eq. 9
    tau=(sig1-sig3)*np.sqrt(ds1s3)/(ds1s3+1)
#    4. Generating spline defined by sign vs tau
    sign_tau=UnivariateSpline(sign,tau,s=0)
#    5. returns spline, and minimum and maximum normal stress values
    return sign_tau, sig_t, 5*sig3max


    
def MCstrength_from_HBenv(sign, sign_tau, sign_min, sign_max):
# =============================================================================
#      Input:    
#      sign       1-d array of int or float; values of normal stress range (MPa)
#      sign_tau   UnivariateSpline defined by normal and shear stresses in 
#                      HoekBrownFE_spline()
#      sign_min   float, minimum normal stress used to define sign_tau (MPa)
#      sign_max   float, maximum normal stress used to define sign_tau (MPa)
#          
#      Output:    
#      phi        float; friction angle, (degrees)
#      coh        float; cohesion, (MPa)
#        
#      Description:        
#      Function returns arrays of friction angle (phi), and cohesion (coh)
#      derived from lines tangent to the Hoek-Brown failure envelope. 
#      
# =============================================================================
#    checking if normal stresses are within range
    sign=sign[np.where((sign>sign_min)&(sign<=sign_max))]
#    computing phi and coh
    phi=np.degrees(np.arctan(sign_tau(sign,1)))
    coh=sign_tau(sign)-sign_tau(sign,1)*sign
    return phi,coh











    
