'''
Created on Mar 21, 2018

@author: Keisuke
''
Here I add data for O2
'''

from pylab import *
from McCarty02 import *
from McCarty01 import *
from Savefig import *
import time

def tricho(DiffusionFactor):

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Time setting 
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Tmax=86400 #(s) maximum time step (177-32)
    dT=1800    #(s) time step    (177-32)
    dT=50
    T=arange(0.,Tmax+10**(-10),dT)  #(s) time array (177-32)
    U=T/dT                   #(s) U array for for loop (from a704_16)
    Th=T/3600       #(h) time 
    LD=ones(size(T))  #Light or dark: 1: light, 0: dark
    dl=43200            #(hours) Day length
    LD[T>=dl]=0
    LD[0]=0
    #-----------------------------------
    #photosynthesis related parameters
    #-----------------------------------
    I=400*ones(size(T)) #*sin(pi*T/dl)
    #I=200*sin(pi*T/dl)
    I[LD==0]=0
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll (from a 704)
    Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990) (from a704-16 (in Healey 1985 model))
    Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion) (from a704-16)
    Pmax = 0.002681321
    O=0.01  #(m2 / umol photon) absorption cross section (from a704_16)
    Tau=1    #(s) handling time (from a704_16) (Irina: half life of RNA is less than 1 minute (30 (s)) or so)
    PI=Pmax*(1-exp(-O*I*Tau))     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
    Rd=0.2*ones(size(T)) #* sin(pi*T/dl)     #(dimensionless) The ratio of diazocyte
    Rd[(T>10800) & (T<32400)]=0.45
    Rp=1-Rd     #(dimensionless) The ratio of photosynthetic cells
    Chlmaxg=0.048     #(gchl gC-1) Chlorophyll max of the cell(186-40)
    Chlmaxmol=Chlmaxg*12*55/868      #(molCchl molC-1) Chlorophyll max of the cell in mol (186-40)
    ro=18333    #(molC / m3) carbon density in the cell
    Chl=Chlmaxmol*0.8*ro        #(molCchl m-3) Chlorophyll concentration
    Chl=Chl*Rp
    Ync=1/5
    Ync=1/6.3   #(molN / molC) ratio of N:C  
    Mumax=0.25/86400  #(s-1) maximum growth rate
    lmax=ro*Mumax  #(molC m-3 s-1)
    Kc=ro/5     #(molC / m3) half saturation constant of carbon storage use
    Kn=Kc*Ync   #(molN / m3) half saturation constant of nitrogen storage use
    ResMax=lmax*200  #(molC m-3 s-1) Maximum respiration rate  
    NfixFull=lmax*Ync*25       #(molN m-3 s-1) Maximmum nitrogen fixation rate
    Fnitroge = 0.4
    NfixMax=NfixFull*Rd*Fnitroge
    NstoMax=ro*Ync      #(molN m-3) Maximum nitrogen storage
    CstoMax=ro*2
    NtoC=McCarty02(0.6)
    E = McCarty01(0.6)
    FvFmP=0.5           #(dimensionless) Fv:Fm of photosynthetic cells
    FvFmD=0.1           #(dimensionless) Fv:Fm of nitrogen fixing cells
    FvFm=FvFmP*Rp + FvFmD*Rd
    
    #Output for paper parameter values==============
#     print(ResMax)
#     print(NfixFull)
#     print(NstoMax)
#     print(lmax)
#     print(O*Tau)
#     print(Pmax)
#     print(Kc)
#     print(CstoMax)
    #---------------------------------------
    #Array preparation for each prarameter
    #---------------------------------------
    o=copy(T)*0             # this creates zero array for the right size for the time steps
    dCsto=copy(o)
    dNsto=copy(o)
    dOx1=copy(o)
    dOx2=copy(o)
    dOx3=copy(o)
    Pchl=copy(o)
    PI=copy(o)
    Res=copy(o)
    ls=copy(o)
    Csto=copy(o)
    Nfix=copy(o)
    Res1=copy(o)
    Res2=copy(o)
    Nsto=copy(o)
    Ox1=copy(o)
    Ox2=copy(o)
    Res1=copy(o)
    Cn2fix=copy(o)
    Fox13=copy(o)
    Fox12=copy(o)
    Fox32=copy(o)
    Fox34=copy(o)
    Ox1=copy(o)
    Ox2=copy(o)
    Ox3=copy(o)
    Csto2=copy(o)
    CcN2fix=copy(o)
    Res2ex=copy(o)
    ls1=copy(o)
    ls2=copy(o)
    Res2c=copy(o)
    which=copy(o)
    Nf1=copy(o)
    Nf2=copy(o)
    Nf3=copy(o)
    ResN2fix=copy(o)
    ResPro=copy(o)
    
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #Initial condition
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    Ox1[0]=0.213
    Ox2[0]=0.213
    Ox3[0]=0.213
    Csto[0]=1000
    Nsto[0]=1000
        
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For statement
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    for i in U:
        Ca0=1.8
        Ca=5751/2    #Constant for A
        #Ca=100
        #E=0.7
        
        A13=2     *Ca * DiffusionFactor
        A12=0.0345    *Ca * DiffusionFactor
        A32=2    *Ca * DiffusionFactor
        A34=2      *Ca0 
        
        A13=A13*Rp[i]
        A32=A32*Rd[i]
        
        Ox4=0.213
        OxCri=0.1
        Yon=8/4.28      #(mol O2/ molN2) nitrogen fixation to oxygen conversion factor
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #Preparing parts
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        PI[i] = Pmax * (1 - exp(-O * I[i] * Tau)) * (CstoMax-Csto[i])/CstoMax     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
        
        Clim=Csto[i] / (Csto[i] + Kc)
        Nlim=Nsto[i] / (Nsto[i] + Kn)
        ls[i] = lmax * min(Clim, Nlim)    #(C mol m-3 s-1) biomass synthesis rate
        
        if Clim>Nlim:
            which[i]=1
        else:
            which[i]=0
        
        Res1[i] = ls[i] * E
        Csto2[i] = Csto[i] / Rd[i]
        Res2c[i] = ResMax*Csto2[i] / (Csto2[i] + Kc) * Rd[i]
        
        #=========================================================================================================
        # O2 flux part Kei210-47 + 02EquationSolving02simplified.nb and 03EquationSolving03different.nb
        #=========================================================================================================
        
        #-------------------------------------------------------
        # When Ox2 = 0 "03EquationSolving03different.nb"
        #-------------------------------------------------------
        a = -Res1[i]  + PI[i] * Chl[i]
        Ox1[i] = -((-a*A13-a*A32-a*A34-A13*A34*Ox4)/(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34))
        Ox3[i] = -((-a*A13-A12*A34*Ox4-A13*A34*Ox4)/(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34))
        b   = -(a*A12*A13+a*A12*A32+a*A13*A32+a*A12*A34+A12*A13*A34*Ox4+A12*A32*A34*Ox4+A13*A32*A34*Ox4)\
               /(A12*A13+A12*A32+A13*A32+A12*A34+A13*A34)
        
        Ox2[i] = 0
        
        Res2[i] = -b
        
        #-------------------------------------------------------
        # Testing b
        #-------------------------------------------------------
        
        if Res2[i] > Res2c[i]:
        
            #-------------------------------------------------------
            # When Ox2 =/ 0 "02EquationSolving02simplified.nb"
            #-------------------------------------------------------
            b= -Res2c[i]
            Ox1[i] = -((-a*A12*A13-a*A12*A32-a*A13*A32-a*A12*A34-a*A32*A34-A12*A13*b-A12*A32*b-A13*A32*b\
                     -A12*A34*b-A12*A13*A34*Ox4-A12*A32*A34*Ox4-A13*A32*A34*Ox4)/((A12*A13+A12*A32+A13*A32)*A34))
            Ox2[i] = -((-a*A12*A13-a*A12*A32-a*A13*A32-a*A12*A34-A12*A13*b-A12*A32*b-A13*A32*b-A12*A34*b\
                     -A13*A34*b-A12*A13*A34*Ox4-A12*A32*A34*Ox4-A13*A32*A34*Ox4)/((A12*A13+A12*A32+A13*A32)*A34))
            Ox3[i] = (a+b+A34*Ox4)/A34
        
            Res2[i] = Res2c[i]
        
        #=========================================================================================================
        
        Nf1[i] = max((OxCri - Ox2[i]) / OxCri, 0)
        Nf2[i] = max((NstoMax - Nsto[i]) / NstoMax, 0)
        Nf3[i] = Csto2[i] / (Csto2[i] + Kc)
        
        Nfix[i] = NfixMax[i] * Nf1[i] * Nf2[i] * Nf3[i]   #(molN m-3 s-1) Nitrogen fixation rate per biomass carbon
        Nfix[LD==0] = 0
        CcN2fix[i] = Nfix[i]         #(molC m-3 s-1) carbon consumption for nitrogen fixation 
        ResN2fix[i] = Nfix[i] * NtoC            #(molC m-3 s-1) Respiration for nitrogen fixation (in carbon)
        ResPro[i] = Res2[i] - ResN2fix[i] 
        Res[i] = Res1[i] + Res2[i]
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #time step change
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        dCsto[i] = PI[i] * Chl[i] - Res[i] - ls[i] - Cn2fix[i] - CcN2fix[i]
        dNsto[i] = Nfix[i] - ls[i] * Ync
        
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        #obtaining next time step values
        #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        if i<U[-1]:
            Csto[i+1] = Csto[i] + dCsto[i] * dT
            Nsto[i+1] = Nsto[i] + dNsto[i] * dT
        
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        
        UnitConversion = 1000
    return mean(Ox1[:dl/dT])*UnitConversion, mean(Ox2[:dl/dT])*UnitConversion, mean(Ox3[:dl/dT])*UnitConversion  
    
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

rcParams.update({'lines.markersize': 10})
rcParams.update({'lines.markeredgewidth': 0.5})
rcParams.update({'font.size': 24})
rcParams.update({'lines.linewidth': 3.5})
rcParams.update({'figure.autolayout': True})
rcParams['figure.figsize']=8,6.5
rcParams.update({'figure.facecolor':'W'})
     
#   rcParams.update({'mathtext.default': 'regular' })
rcParams.update({'patch.edgecolor':'none'})
rcParams.update({'xtick.major.pad': 10})
rcParams.update({'ytick.major.pad': 10})

rcParams.update({'axes.linewidth':2})
rcParams.update({'xtick.major.width':1.5})
rcParams.update({'ytick.major.width':1.5})
#OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

Dstep = 0.01
Dmin = -4
Dmax = 0
DiffusionFactorArray = arange(Dmin,Dmax + Dstep,Dstep)
DiffusionFactorArray = 10**DiffusionFactorArray
Dmin = 10**Dmin
Dmax = 10**Dmax

print(DiffusionFactorArray)
O = zeros(size(DiffusionFactorArray))
Ox1array = copy(O)
Ox2array = copy(O)
Ox3array = copy(O)
Ox4array = copy(O) + 213


i = 0
t00 = time.time()
t =time.time()
for DiffusionFactor in DiffusionFactorArray:    
    t0 = t
    Ox1array[i],Ox2array[i],Ox3array[i] = tricho(DiffusionFactor)
    t=time.time()
    print(i,round(t-t0,2))
    i = i + 1

UnitConversion = 1000

figure(1)
semilogx(DiffusionFactorArray,Ox1array,label='[O$_2$]$_1$')
semilogx(DiffusionFactorArray,Ox2array,label='[O$_2$]$_2$')
semilogx(DiffusionFactorArray,Ox3array,label='[O$_2$]$_3$')
semilogx(DiffusionFactorArray,Ox4array,label='[O$_2$]$_4$')
semilogx(DiffusionFactorArray,303*ones(size(DiffusionFactorArray)),':',color='r',label='[O$_2$]$_{3 data}$')

# #=======================
# # Error shading
# #=======================
# Error = 51.55434
# y1 = 318*ones(size(DiffusionFactorArray)) + Error
# y2 = 318*ones(size(DiffusionFactorArray)) - Error
# fill_between(DiffusionFactorArray, 318*ones(size(DiffusionFactorArray)) + Error, 318*ones(size(DiffusionFactorArray)) - Error, where=y1 >= y2, facecolor='#FDD2D0',edgecolor = "none")
# #=======================

xlabel('Relative diffusivity',labelpad=8)
ylabel('O$\mathregular{_2}$ (uM)',labelpad=8)
#legend(loc='upper right',borderaxespad=0.5,ncol=2, fontsize=20)
xlim(Dmin - 1e-17,Dmax + 1e-17 )
ylim(0,800)
ylim(ymin=-600/35)

#Savefig######################
First='Tricho\\'
RunName='004504'
Folderloc=First+RunName
Savefig(Folderloc,'O2',600)

##############################

figure(2)
semilogx(DiffusionFactorArray,Ox1array,label='[O$_2$]$_1$')
semilogx(DiffusionFactorArray,Ox2array,label='[O$_2$]$_2$')
semilogx(DiffusionFactorArray,Ox3array,label='[O$_2$]$_3$')
semilogx(DiffusionFactorArray,Ox4array,label='[O$_2$]$_4$')
xlabel('Relative diffusivity',labelpad=8)
ylabel('O$\mathregular{_2}$ (uM)',labelpad=8)
#legend(loc='upper right',borderaxespad=0.5, fontsize=23)
xlim(Dmin - 1e-17,Dmax + 1e-17 )

#Savefig######################
Folderloc=First+RunName
ylim(ymin=-600/35)
Savefig(Folderloc,'O2-NoLegend',600)

##############################

print('Total time', round(time.time()-t00,2),'(s)')
show()


    