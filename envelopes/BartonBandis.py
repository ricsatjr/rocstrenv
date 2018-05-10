# =============================================================================
# DESCRIPTION:
# Module computing shear strength envelope from a given set of 
# Barton-Bandis strength parameters and a range of normal stresses.   
# 
# REFERENCE
# Barton, N. R., & Bandis, S. C. (1990). Review of predictive capabilities 
# of JRC-JCS model in engineering practice. In Rock Joints, Proc int symp on 
# rock joints, Loen, Norway (eds N. Barton and O. Stephenson) (pp. 603-610).
# =============================================================================


# libraries and modules
import numpy as np

def BBenv(sig_n,JRC0,JCS0,phi_r,Ln):
# =============================================================================
#     Input:    
#     sig_n     array of floats; normal stresses (MPa)
#     JRC0       int; Joint roughness coefficient (unitless), measured at 10cm scale
#     JCS0       int or float; uniaxial compressive strength of joint (MPa), measured at 10cm scale
#     phi_r     float; residual friction angle of joint (degrees)
#     Ln        float; length of joint being assessed (m)
#     
#     Output:    
#     tau    array of floats; shear strength corresponding to sig_n

#       
#     Description:        
#     Function returns an array of shear strength values corresponding to  
#     normal stresses, sig_n, as defined by the Barton-Bandis failure criterion. 
#     
#     Variables and equations are from Barton and Bandis (1990).
# =============================================================================
#    removing negative normal stresses
    sig_n=np.where(sig_n>0,sig_n,np.nan)
#    reducing JRC and JCS based on length of joint
    JRCn,JCSn=BBscaling(JRC0,JCS0,Ln)
#    eq. 4 
    Phi=JRCn*np.log10(JCSn/sig_n)+phi_r
    tau=sig_n*np.tan(np.radians(Phi))
    return np.where((Phi<85)&(tau>0),tau,np.nan)

def BBscaling(JRC0,JCS0,Ln):
# =============================================================================
#     Input:    
#     JRC0       int; Joint roughness coefficient (unitless), measured at 10cm scale
#     JCS0       int or float; uniaxial compressive strength of joint (MPa), measured at 10cm scale
#     Ln        float; length of joint being assessed (m)
#     
#     Output:    
#     JRCn       int or float; scaled JRC0 value (unitless)
#     JCSn       int or float; scaled JCS0 value (MPa)
#       
#     Description:        
#     Function returns scaled values of JRC0 and JCS0 as a function of joint length, Ln
#     
#     Variables and equations are from Barton and Bandis (1990).
# =============================================================================
    L0=0.1 #m
#    eq. 9
    JRCn=JRC0*(Ln/L0)**(-0.02*JRC0)
#    eq. 10
    JCSn=JCS0*(Ln/L0)**(-0.03*JRC0)
    return JRCn,JCSn



    



