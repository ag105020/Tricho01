'''
Created on Mar 21, 2018

@author: Keisuke
''
Beginning of the modeling
Model largely based on Kei 210-27~
'''

from pylab import *
from McCarty00 import *
from Savefig import *

def tricho():

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
    I=200*ones(size(T)) #*sin(pi*T/dl)
    I[LD==0]=0
    Mchl=893.49             #(g / mol chlorophyll) mollar mass of chlorophyll (from a 704)
    Pmax0=7                     #(g C /(g Chl h) Maximum production rate per chlorophyll (around 6 by Cullen 1990) (from a704-16 (in Healey 1985 model))
    Pmax=Pmax0*Mchl/12/3600/55       #(mol C s-1 mol chl-1) carbon fixing rate (156-10) (156-15) for unit conversion) (from a704-16)
    O=0.01  #(m2 / umol photon) absorption cross section (from a704_16)
    Tau=5    #(s) handling time (from a704_16) (Irina: half life of RNA is less than 1 minute (30 (s)) or so)
    PI=Pmax*(1-exp(-O*I*Tau))     #(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)
    Rd=0.2*ones(size(T)) #* sin(pi*T/dl)     #(dimensionless) The ratio of diazocyte
    Rd[(T>10800) & (T<32400)]=0.4
    Rp=1-Rd     #(dimensionless) The ratio of photosynthetic cells
    Chlmaxg=0.048     #(gchl gC-1) Chlorophyll max of the cell(186-40)
    Chlmaxmol=Chlmaxg*12*55/868      #(molCchl molC-1) Chlorophyll max of the cell in mol (186-40)
    ro=18333    #(molC / m3) carbon density in the cell
    Chl=Chlmaxmol*0.8*ro        #(molCchl m-3) Chlorophyll concentration
    Chl=Chl*Rp
    Ync=1/5
    Ync=1/5.5   #(molN / molC) ratio of N:C  
    Mumax=0.2/86400  #(s-1) maximum growth rate
    lmax=ro*Mumax  #(molC m-3 s-1)
    Kc=ro/5     #(molC / m3) half saturation constant of carbon storage use
    Kn=Kc*Ync   #(molN / m3) half saturation constant of nitrogen storage use
    ResMax=lmax*200  #(molC m-3 s-1) Maximum respiration rate  
    NfixFull=lmax*Ync*10       #(molN m-3 s-1) Maximmum nitrogen fixation rate
    NfixMax=NfixFull*Rd
    NstoMax=ro*Ync      #(molN m-3) Maximum nitrogen storage
    CstoMax=ro*2
    NtoC=McCarty00()
    FvFmP=0.5           #(dimensionless) Fv:Fm of photosynthetic cells
    FvFmD=0.1           #(dimensionless) Fv:Fm of nitrogen fixing cells
    FvFm=FvFmP*Rp + FvFmD*Rd
    
    #Output for paper parameter values==============
    print(ResMax)
    print(NfixFull)
    print(NstoMax)
    print(lmax)
    print(O*Tau)
    print(Pmax)
    print(Kc)
    print(CstoMax)
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
    Ox1[0]=0.225
    Ox2[0]=0.225
    Ox3[0]=0.225
    Csto[0]=1000
    Nsto[0]=1000
        
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    #For statement
    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    for i in U:
        
        Ca=2
        #Ca=100
        E=0.7
        
        A13=2     *Ca
        A12=0.2    *Ca
        A32=2    *Ca
        A34=2      *Ca
        
        A13=A13*Rp[i]
        A32=A32*Rd[i]
        
        Ox4=0.225
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
        # O2 flux part Kei210-33 + 02EquationSolving02simplified.nb and 03EquationSolving03different.nb
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
    
    Dh=ls/ro*3600   #(h-1) growth rate 
    Ph=PI*Chl/ro*3600   #(h-1) photosynthesis rate per carbon
    Rh=Res/ro*3600   #(h-1) respiration rate per carbon
    
    Dd=ls/ro*86400    #(d-1) growth rate
    Pd=PI*Chl/ro*86400    #(d-1) photosynthesis rate per carbon
    Rday=Res/ro*86400    #(d-1) respiration rate per carbon
    R1d=Res1/ro*86400  #(d-1) respiration rate in cell1
    R2d=Res2/ro*86400  #(d-1) respiration rate in cell2
    
    NfdC=CcN2fix/ro*86400 #(molC molC-1 d-1) nitrogen fixation rate 
    Nfd=Nfix/ro*86400 #(molN molC-1 d-1) nitrogen fixation rate
    
    dNd=dNsto/ro*86400 #(d-1) change in nitrogen storge
    
    RNFd=ResN2fix/ro*86400  #(d-1) respiration for nitrogen fixation
    RPd = ResPro/ro*86400 #(d-1) respiratory protection
    CSTOd=dCsto/ro*86400    #(d-1) storage acumulation rate
    CSTOdP=copy(CSTOd); CSTOdP[CSTOd<=0]=0    #(d-1) positive part of cSTOd
    CSTOdN=copy(CSTOd); CSTOdN[CSTOd>=0]=0
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #Plot
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    rcParams.update({'lines.markersize': 10})
    rcParams.update({'lines.markeredgewidth': 0.5})
    rcParams.update({'font.size': 22})
    rcParams.update({'lines.linewidth': 3.5})
    rcParams.update({'figure.autolayout': True})
    rcParams['figure.figsize']=8,6.5
    rcParams.update({'figure.facecolor':'W'})
         
    rcParams.update({'mathtext.default': 'regular' })
    rcParams.update({'patch.edgecolor':'none'})
    rcParams.update({'xtick.major.pad': 8})
    rcParams.update({'ytick.major.pad': 8})

    #OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

    fig1=1
    fig2=1
    fig3=1
#     fig4=1
#     fig5=1
#     fig6=1
#     fig7=1
#     fig8=1
#     fig9=1
#     fig10=1
    
    paperplot=0
    
    Xt='Time (h)'
    Xmax=Tmax/3600
    Xstep=3
    
    First='Tricho\\'
    RunName='002800'
    Folderloc=First+RunName
    
    if paperplot!=1:
        if fig1==1:
            figure(1,figsize=(13,10))
            subplot(2,2,1)
            plot(Th,Rd)
            title('Rd')
            ylim(0,1)
            ylabel('Rd')
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            
            subplot(2,2,2)
            plot(Th,Dd)
            title('Growth rate (d$^{-1}$)')
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            ylabel('Growth rate (d$^{-1}$)')
            
            subplot(2,2,3)
            plot(Th,Pd,label='Photo')
            plot(Th,NfdC,label='Nfix')
            plot(Th,Rday,label='Resp')
            title('C Flux')
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            ylabel('C Flux (mol C mol C$^{-1}$ d$^{-1}$)')
            legend(loc=1, fontsize=20)
            
            subplot(2,2,4)
            plot(Th,Csto,label='Csto')
            plot(Th,Nsto,label='Nsto')
            title('Storage (mol m$^{-3}$)')
            legend(loc=1,fontsize=20)
            ylabel('Sto (mol m$^{-3}$)')
            ylim(ymin=-1000)
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            
            Pad=2; W_pad=2; H_pad=2
            tight_layout(pad=Pad,w_pad=W_pad,h_pad=H_pad)
        
        if fig2==1:
            figure(2,figsize=(13,10))
            subplot(2,2,1)
            plot(Th,Rday,label='Rday')
            plot(Th,R1d,label='Rd1')
            plot(Th,R2d,label='Rd2')
            plot(Th,Pd,label='Pd')
            legend()
            
            subplot(2,2,2)
            plot(Th,Pd)
            title('Pd')
            
            subplot(2,2,3)
            plot(Th,dNd)
            title('dNd')
            
            subplot(2,2,4)
            plot(Th,which)
            title('0=Clim, 1=Nlim')
            ylim(-0.2,1.2)
            yticks([0,1])
        
        if fig3==1:
            figure(3,figsize=(13,10))
            subplot(2,2,1)
            plot(Th,Nfd)
            title('Nfix')
            ylabel('Nfix (mol N mol C$^{-1}$ d$^{-1}$)')
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            
            ax=subplot(2,2,2)
            plot(Th,Ox1,label='[O$_2$]$_1$')
            plot(Th,Ox2,label='[O$_2$]$_2$')
            plot(Th,Ox3,label='[O$_2$]$_3$')
            plot(Th,Ox4*ones(size(Th)),label='[O$_2$]$_4$')
            legend(loc=1,borderaxespad=0.5, fontsize=20)
           
            title('O$_2$ concentrations', y=1.02)
            ylim(ymin=-0.01)
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            ylabel('O$_2$ (mM)')
            
            subplot(2,2,3)
            Names = ["Bp","ResBp","Nfix","ResN2","RP","Csto"]
            Colors = ['blue','orange','cyan','red','green','greenyellow']
            for i in arange(size(Names)):
                plot([],[],linewidth=10,color=Colors[-i-1],label=Names[-i-1])
            stackplot(Th,Dd,R1d,Nfd,RNFd,RPd,CSTOdP,colors=Colors)
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            xlim(xmin=0)
            ylabel('C flux (mol C mol C$^{-1}$ d$^{-1}$)')
            title('Fate of C')
            legend(fontsize=20)
            
            subplot(2,2,4)
            plot(Th,FvFm)
            title('Fv/Fm')
            xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
            ylabel('Fv/Fm')
            ylim(0,0.600001)
            Pad=2; W_pad=2; H_pad=2
            tight_layout(pad=Pad,w_pad=W_pad,h_pad=H_pad)
         
        show(block=True)
    if paperplot==1:
        
        #########################
        rcParams.update({'font.size': 28})
        ######################### 
        figure(5, figsize=(16,13))
        subplot2grid((2,4),(0,0),colspan=2)
        plot(Th,Pd,label='Photo')
        plot(Th,NfdC,label='Nfix')
        plot(Th,Rday,label='Resp')
        title('C Flux')
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        ylabel('C Flux (mol C mol C$^{-1}$ d$^{-1}$)')
        legend(loc=1, fontsize=28)    
        
        subplot2grid((2,4),(0,2),colspan=2)
        plot(Th,Nfd)
        title('Nfix')
        ylabel('Nfix (mol N mol C$^{-1}$ d$^{-1}$)')
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        
        subplot2grid((2,4),(1,1),colspan=2)
        plot(Th,FvFm)
        title('Fv/Fm')
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        ylabel('Fv/Fm')
        ylim(0,0.600001)
        Pad=2; W_pad=2; H_pad=2
        tight_layout(pad=Pad,w_pad=W_pad,h_pad=H_pad)
        
        Savefig(Folderloc,5,300)
        
        
        figure(6)
        plot(Th,Dd)
        title('Growth rate (d$^{-1}$)')
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        ylabel('Growth rate (d$^{-1}$)')
        Savefig(Folderloc,6,300)
        
        figure(7)
        plot(Th,Csto,label='Csto')
        plot(Th,Nsto,label='Nsto')
        title('Storage (mol m$^{-3}$)')
        legend(loc=1,fontsize=22)
        ylabel('Sto (mol m$^{-3}$)')
        ylim(ymin=-1000)
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        Savefig(Folderloc,7,300)
        
        figure(8)
        Names = ["Bp","ResBp","Nfix","ResN2","RP","Csto"]
        Colors = ['blue','orange','cyan','red','green','greenyellow']
        for i in arange(size(Names)):
            plot([],[],linewidth=15,color=Colors[-i-1],label=Names[-i-1])
        stackplot(Th,Dd,R1d,Nfd,RNFd,RPd,CSTOdP,colors=Colors)
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        xlim(xmin=0)
        ylabel('C flux (mol C mol C$^{-1}$ d$^{-1}$)')
        title('Fate of C')
        legend(fontsize=25)
        Savefig(Folderloc,8,300)
        
        figure(9)
        SC=Dd+R1d+Nfd+RNFd+RPd+CSTOdP  #sum of carbon
        Names = ["Bp","ResBp","Nfix","ResN2","RP","Csto"]
        Colors = ['blue','orange','cyan','red','green','greenyellow']
        for i in arange(size(Names)):
            plot([],[],linewidth=15,color=Colors[-i-1],label=Names[-i-1])
        stackplot(Th,Dd/SC,R1d/SC,Nfd/SC,RNFd/SC,RPd/SC,CSTOdP/SC,colors=Colors)
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        xlim(xmin=0)
        ylim(ymax=1.0000001)
        ylabel('C flux ($\%$)')
        title('Fate of C (in $\%$)')
        #legend(fontsize=25)
        Savefig(Folderloc,9,300)
        
        figure(10)
        plot(Th,Ox1,label='[O$_2$]$_1$')
        plot(Th,Ox2,label='[O$_2$]$_2$')
        plot(Th,Ox3,label='[O$_2$]$_3$')
        plot(Th,Ox4*ones(size(Th)),label='[O$_2$]$_4$')
    
        title('O$_2$ concentrations', y=1.04)
        ylim(ymin=-0.01)
        xlabel(Xt); xlim(0,Xmax); xticks(arange(0,Xmax+Xstep,Xstep))
        ylabel('O$_2$ (mM)')
        ylim(ymax=0.7000001)
        if Ca>50:
            Savefig(Folderloc,'10 High Diffusion',300)
        else:
            Savefig(Folderloc,10,300)
            legend(loc=1,borderaxespad=0.5, fontsize=25)
        #rcParams.update({'font.size': 28})
    
        show()
    
    
tricho()
    
    
    
    