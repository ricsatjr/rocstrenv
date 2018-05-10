import numpy as np
import matplotlib.pyplot as plt
import HoekBrown as HB
import BartonBandis as BB

#MAIN

plt.close('all')

H=50.   #slope height, meters
fig,ax=plt.subplots(nrows=1,ncols=2,sharex=True,sharey=True)


label1=['a','b','c']
GSI=np.array([15,30,60])
mi=np.array([5,15,25])
sig_ci=np.array([5.,20.,35.])   #MPa
uw=np.array([24.,24.,24.])/1000.  #MPa
D=0.7




for o in range(len(sig_ci)):
    sign_tau, sig_t,sign_max=HB.HBenv(GSI[o],D,mi[o],sig_ci[o],uw[o],H)
    sig_n=np.linspace(sig_t,0.5*uw[o]*9.81*H,1000)
    tau=sign_tau(sig_n)
    ax[0].set_title('Hoek-Brown')
    ax[0].plot(sig_n,tau, label='GSI='+str(GSI[o])+' mi='+str(mi[o])+' UCS='+str(sig_ci[o]))
    ax[0].set_xlabel('normal stress, MPa')
    ax[0].set_ylabel('shear stress, MPa')
    ax[0].grid(b=None, which='both', axis='both',ls=':')
    ax[0].legend(fontsize='small',loc='upper left')

##############################################################
    
label2=['a','b','c']
JRC=[16.9,0.5,8.9]
res_phi=[29.,26.,28.]   #deg
JCS=[96.,50.,92.]        #MPa
uw=np.array([24.,24.,24.])/1000.  #MPa
Ln=0.1                  #meters

for o in range(len(JCS)):
    sig_n=np.linspace(0,0.5*uw[o]*9.81*H,1000)
    tau=BB.BBenv(sig_n,JRC[o],JCS[o],res_phi[o],Ln)
    ax[1].set_title('Barton-Bandis')
    ax[1].plot(sig_n,tau, label='JRC='+str(JRC[o])+' phi_r='+str(res_phi[o])+' JCS='+str(JCS[o]))
    ax[1].set_xlabel('normal stress, MPa')
    ax[1].grid(b=None, which='both', axis='both',ls=':')
    ax[1].legend(fontsize='small',loc='upper left')
plt.show()
