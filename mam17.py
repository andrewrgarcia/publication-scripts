# -*- coding: utf-8 -*-
#mam17.py
#MAM article figures

from mamdata import *
from scipy.stats import tstd
from decimal import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sigfig
import scipy.stats as spstats
import pylab 

def fig():
## INDEX ##
    print (70 * '-')
    print ("          Script for MAM Article Figures")
    print ("          Andrew Garcia, 2017")
    print (70 * '.')
    print ("Figure ID        Figure Description ")
    print (70 * '-')
    print ("f1:              Microsphere size histograms from microsphere preparations with varying ")
    print ("                 shear rates made with alginate volume fractions of 0.984% v/v and using")
    print ("                 a high viscosity continuous phase")
    print ("f2:              Microsphere size histograms from microsphere preparations with varying ")
    print ("                 alginate volume fractions made under applied shear rates of 9400 s-1 ")
    print ("                 and using a high viscosity continuous phase")
    print ("f3:              Microsphere size reduction with increasing applied shear rates observed ")
    print ("                 from magnetic alginate microsphere preparations made with alginate ")
    print ("                 volume fractions of 0.984% v/v and using a high viscosity continuous phase.")
    print ("f4:              Microsphere size increase with increasing alginate volume fractions on ")
    print ("                 microsphere size observed from magnetic alginate microsphere preparations ")
    print ("                 made under applied shear rates of 9400 s-1 and using a high viscosity ")
    print ("                 continuous phase.")
    print ("f5:              Microsphere size histograms from microsphere preparations with varying ")
    print ("                 shear rates made with alginate volume fractions of 0.984% v/v and using")
    print ("                 a low viscosity continuous phase")
    print ("f6:             Microsphere size histograms from microsphere preparations with different ")
    print ("                 magnetic nanoparticle concentrations, made under applied shear rates of") 
    print ("                 9400 s-1, with alginate volume fractions of 0.984% v/v and using a high")
    print ("                 viscosity continuous phase.")
    print (70 * '-')
        
    study_handle = raw_input('Enter Figure ID:')
    study = study_handle
    
    '''Define Global Variables:'''
    
    global shear 
    [micsM, meanV, sigmaV, micsMoct, meanVoct, sigmaVoct,dispM,\
    meanVdisp, sigmaVdisp, bangsM, meanVangs, sigmaVangs, mudxM, meanVmudx, sigmaVmudx]=data()
    
    d32moil=np.zeros(8)
    for i in range(0,8):
        d32moil[i]=sauter(micsM[i])
        
    d32oct=np.zeros(8)
    for i in range(0,8):
        d32oct[i]=sauter(micsMoct[i])
        
    d32mudx=np.zeros(7)
    for i in range(0,7):
        d32mudx[i]=sauter(mudxM[i])
        
    
    if study=='f2':
        
        mu_N=np.array([ 100.,  100.,  101.,  100.,  100.,  100.,  100.])
        
        '''HISTOGRAMS'''
        bins =20
        # A= SUM(a_ij)
        histM = [[0 for j in range(7)] for i in range(np.shape(mudxM[2])[0])]
        binsM = [[0 for j in range(7)] for i in range(np.shape(mudxM[2])[0])]
        # Create histograms and normalize total count to 1 (data):
        for i in range(0,7):
            histM[i],binsM[i]= np.histogram(mudxM[i], bins = bins)
            histM[i] = [ float(n)/sum(histM[i]) for n in histM[i]]
        # Set histograms' parameters (data)
        center = [[0 for j in range(7)] for i in range(20)]
        width=np.zeros(7)
        for i in range(0,7):
            center[i] = (binsM[i][:-1]+binsM[i][1:])/2
            width[i]= 1*(binsM[i][1]-binsM[i][0]) 
            
        f, ax = plt.subplots(7)            
        # Generate histograms + hist. annotations/description + hist. color style  
        
        'normality tests (p-value)'
        Np=np.zeros(7)
        for i in range(0,7):
            Np[i]=spstats.mstats.normaltest(mudxM[i])[1]
            
        #normality tests
        #for k in range(0,7):
        #    #ax[k].text(10, 0.22,r"$p$="+str(round(Np[k],5)))
        #    ax[k].text(10, 0.22,r"$p$="+str("{:.2e}".format(Decimal(Np[k]))))
                    
        cmap = mpl.cm.cool
        for k in range(0,7):
            #ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color=cmap((2*k) / float(11)))
            ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color='w')
            #ax[k].text(10,0.22,"$d_{32}$="+str(sigfig.round_sig(d32mudx[k], 2))+r"$\; \mathrm{\mu m}$")
            ax[k].set_ylabel(''), ax[k].set_xscale('log'), ax[k].set_xlim([0.1, 100]),ax[k].set_ylim([0, 0.3])
        
        #for k in range(0,7):
            #ax[k].text(10,0.22,"$d_{32}$="+str(sigfig.round_sig(d32mudx[k], 2))+r"$\; \mathrm{\mu m}$")
        
        '''non-uniform distribution - measures of dispersion'''
        'median'
        medv=np.zeros(7)
        for i in range(0,7):
            medv[i]=np.median(mudxM[i])
    
        'interquartile range'
        Q1=np.zeros(7)
        for i in range(0,7):
            Q1[i]=IQR(a=mudxM[i])[1]
            
        Q3=np.zeros(7)
        for i in range(0,7):
            Q3[i]=IQR(a=mudxM[i])[2]
        
        #diameter sig figs
        for k in range(0,5):
            #ax[k].text(10, 0.12,r"$d$="+str(round(meanVmudx[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaVmudx[k], 1))+r"$\; \mathrm{\mu m}$")
            #ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1))+r" $IQR$=["+str(round(Q1[k],1))+r","+str(round(Q3[k],1))+"]"+r"$\; \mathrm{\mu m}$")
            ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1)))
        for k in range(5,7):
            #ax[k].text(10, 0.12,r"$d$="+str(round(meanVmudx[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaVmudx[k], 1))+r"$\; \mathrm{\mu m}$")
            #ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1))+r" $IQR$=["+str(round(Q1[k],1))+r","+str(round(Q3[k],1))+"]"+r"$\; \mathrm{\mu m}$")
            ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1)))
    
        #disp phase visc sig figs
        C = np.array([20,30,40,50,60,70,80])
        #Alginate in dispersed phase, Algd?
        #Calg=C/2
        
        mudV=np.zeros(7)
        for i in range(0,7):
            #Ns/m2
            mudV[i]=0.001*einstein(C=C[i])
            #mPa-s
            #mudV[i]=einstein(C=C[i])
        phi = 100*np.array([0.00494364247577615,0.00739717920899497,0.00983864620228257,0.0122681322995387,0.0146857254748385,0.0170915128430511,0.019485580670304])
    
    
        for k in range(0,7):
            ax[k].text(10,0.22,"$\phi$="+str(sigfig.round_sig(phi[k],3))+r" $\%$")
            ax[k].text(0.11,0.22,"$N$="+str(mu_N[k]))
            #ax[k].text(0.11,0.20,"$\mu_d$="+str(sigfig.round_sig(mudV[k],3))+r"$\; \mathrm{mPa-s}$")
            #ax[k].text(0.11,0.20,"$C_{alg,d}$="+str(sigfig.round_sig(Calg[k],2))+r"$\; \mathrm{mg/mL}$")
        
        ax[6].set_xlabel('Microsphere Size, $d$ ($\mathrm{\mu m}$)',size=13)
        
        
        # Fine-tune figure
        f.subplots_adjust(hspace=0,wspace=0,left=0.08,right=0.95)
        plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
        plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
        plt.setp(ax[6].get_xticklabels(), visible=True)
        #plt.suptitle(r"Alginate in Dispersed Phase, $C_{alg,d}$ ($\mathrm{mg/mL}$)",size=13)
        plt.suptitle(r"Alginate Volume Fraction, $\phi$ ($\%$)",size=13)
        #plt.tsubplots_adjust(bottom=0.2)
        plt.show()
    
    #three
    if study=='f1':
        
        '''SHEAR RATES (gamma) vector'''
        RSc=3.82913349489977 
        gamma2=shear(N=16800,Ri=2,R=1.77)
        gamma4=shear(N=8000,Ri=2,R=1.77)
        gamma6=shear(N=1000,Ri=4,R=RSc)
        gamma8=shear(N=300,Ri=4,R=RSc)
        gamma7=shear(N=500,Ri=4,R=RSc)
        gamma5=shear(N=2000,Ri=4,R=RSc)
        gamma3=shear(N=4000,Ri=4,R=RSc)
        gamma1=shear(N=6000,Ri=4,R=RSc)
        gammaV = np.array([gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8])
        gamma_N=np.array([ 103.,  139.,  101.,  114.,  153.,  159.,   82.,  140.])
    
        'REVERSED STACK'
        gammaV=gammaV[::-1]
        gamma_N=gamma_N[::-1]
        micsM=micsM[::-1]
        meanV=meanV[::-1]
        sigmaV=sigmaV[::-1]
        d32moil=d32moil[::-1]
        
        '''HISTOGRAMS'''
        bins =20
        # A= SUM(a_ij)
        histM = [[0 for j in range(8)] for i in range(np.shape(micsM[2])[0])]
        binsM = [[0 for j in range(8)] for i in range(np.shape(micsM[2])[0])]
                
        # Create histograms and normalize total count to 1 (data):
        for i in range(0,8):
            histM[i],binsM[i]= np.histogram(micsM[i], bins = bins)
            histM[i] = [ float(n)/sum(histM[i]) for n in histM[i]]
        # Set histograms' parameters (data)
        center = [[0 for j in range(8)] for i in range(20)]
        width=np.zeros(8)
        for i in range(0,8):
            center[i] = (binsM[i][:-1]+binsM[i][1:])/2
            width[i]= 1*(binsM[i][1]-binsM[i][0])          
                
        'PLOT'
        f, ax = plt.subplots(8)            
        # Generate histograms + hist. annotations/description + hist. color style         
        cmap = mpl.cm.cool
    
        'normality tests (p-value)'
        Np=np.zeros(8)
        for i in range(0,8):
            Np[i]=spstats.mstats.normaltest(micsM[i])[1]
            
        #normality tests
        #for k in range(0,8):
        #    #ax[k].text(10, 0.22,r"$p$="+str(round(Np[k],5)))
        #    ax[k].text(10, 0.22,r"$p$="+str("{:.2e}".format(Decimal(Np[k]))))
        
        for k in range(0,8):
            #ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color=cmap((k+2) / float(11)))
            ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color='w')
            ax[k].set_ylabel(''), ax[k].set_xscale('log'), ax[k].set_xlim([0.1, 100]),ax[k].set_ylim([0, 0.3])
        
        #d32 sauter diameters
        #for k in range(0,8):
            #ax[k].text(10,0.22,"$d_{32}$="+str(sigfig.round_sig(d32moil[k], 2))+r"$\; \mathrm{\mu m}$")
                
        '''non-uniform distribution - measures of dispersion'''
        'median'
        medv=np.zeros(8)
        for i in range(0,8):
            medv[i]=np.median(micsM[i])
    
        'interquartile range'
        Q1=np.zeros(8)
        for i in range(0,8):
            Q1[i]=IQR(a=micsM[i])[1]
            
        Q3=np.zeros(8)
        for i in range(0,8):
            Q3[i]=IQR(a=micsM[i])[2]
        
        'DROP SIZE'
        
        #median nonuniform dist. var.
        for k in range(0,5):
            #ax[k].text(10, 0.12,r"$d$="+str(round(meanVmudx[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaVmudx[k], 1))+r"$\; \mathrm{\mu m}$")
            #ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1))+r" $IQR$=["+str(round(Q1[k],1))+r","+str(round(Q3[k],1))+"]"+r"$\; \mathrm{\mu m}$")
            ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1)))
        for k in range(5,8):
            #ax[k].text(10, 0.12,r"$d$="+str(round(meanVmudx[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaVmudx[k], 1))+r"$\; \mathrm{\mu m}$")
            #ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1))+r" $IQR$=["+str(round(Q1[k],1))+r","+str(round(Q3[k],1))+"]"+r"$\; \mathrm{\mu m}$")
            ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1)))
        
        #diameter sig figs
        #for k in range(0,5):
        #    ax[k].text(10, 0.12,r"$d$="+str(round(meanV[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaV[k], 1))+r"$\; \mathrm{\mu m}$")
        #for k in range(5,8):
        #    ax[k].text(10, 0.12,r"$d$="+str(sigfig.round_sig(meanV[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaV[k], 1))+r"$\; \mathrm{\mu m}$")
        #Shear rate sig figs
        #for k in range(0,2):
        for k in range(0,6):
            #ax[k].text(0.11,0.20,"$S$="+str(sigfig.round_sig(gammaV[k],3))+r"$\; \mathrm{s^{-1}}$")
            ax[k].text(10,0.22,"$\dot{\gamma}$="+str(sigfig.round_sig(gammaV[k],2))+r" $1/s$")
            ax[k].text(0.11,0.22,"$N$="+str(gamma_N[k]))
        #for k in range(2,8):
        for k in range(6,8):
            #ax[k].text(0.11,0.20,"$S$="+str(sigfig.round_sig(gammaV[k],2))+r"$\; \mathrm{s^{-1}}$")
            ax[k].text(10,0.22,"$\dot{\gamma}$="+str(sigfig.round_sig(gammaV[k],3))+r" $1/s$")
            ax[k].text(0.11,0.22,"$N$="+str(gamma_N[k]))
        ax[7].set_xlabel('Microsphere Size, $d$ ($\mathrm{\mu m}$)',size=13)
        
        # Fine-tune figure
        f.subplots_adjust(hspace=0,wspace=0,left=0.08,right=0.95)
        plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
        plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
        plt.setp(ax[7].get_xticklabels(), visible=True)
        plt.suptitle(r"Shear Rate, $\dot{\gamma}$ ($1/s$)",size=13)
        #plt.tsubplots_adjust(bottom=0.2)
        plt.show()
        
        #return micsM
        
        
    if study=='f3':
        
        'shear rates (EXPERIMENTAL)'
        RSc=3.82913349489977 
        gamma2=shear(N=16800,Ri=2,R=1.77)
        gamma4=shear(N=8000,Ri=2,R=1.77)
        gamma6=shear(N=1000,Ri=4,R=RSc)
        gamma8=shear(N=300,Ri=4,R=RSc)
        gamma7=shear(N=500,Ri=4,R=RSc)
        gamma5=shear(N=2000,Ri=4,R=RSc)
        gamma3=shear(N=4000,Ri=4,R=RSc)
        gamma1=shear(N=6000,Ri=4,R=RSc)
        gammaV = np.array([gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8])
    
        
        'LENG AND CALABRESE MODEL'
        
        '''SHEAR RATE'''
        #S=np.linspace(700,14100)
        S=np.linspace(300,22000)
        slengd2=lengTKE(phi=0.9838646,S=S)
                
        """DIAMETER (EXPERIMENTAL DATA) MEANS AND STDEVs"""         
        #mean (shear rate)
        d2=meanV 
        e2=sigmaV      
            
    
        '''MEDIANS'''
        medS=np.zeros(8)
        for i in range(0,8):
            medS[i]=np.median(micsM[i])
            
        boxwidth=0.001*np.ones(7)
                
        
        'SHEAR RATE PLOT'
        
        CIl=np.array([0.25974852526248,0.277306917907307,0.46433653140935,0.707815841779272,0.964512763590039,1.74508220138807,2.96400939983538,4.2795939588573])
        CIu=np.array([1.11648913268801,1.13404752521657,1.3210771377387,1.5645564474285,1.82125336910446,2.60182280919425,3.82075001678743,5.13633459079588])
        
        #CIl=CIl[::-1]
        ##plt.plot(S,slengd2,color='k',label=r'Shear Rate, $\dot{\gamma}$ ($1/s$)')
        plt.plot(S,slengd2,color='k')
        plt.plot(gammaV,medS,marker='o',color='r',linestyle='None')
        'confidence intervals'
        plt.plot(gammaV,CIl,color='magenta')
        plt.plot(gammaV,CIu,color='magenta')
        
        plt.boxplot(micsM,positions=gammaV,widths = 0.001)
        plt.text(4500,20,"$\phi$="+str(0.984)+r' $\%$')
        plt.text(4500,40,"$d = e^{9.22}\phi^{0.75}\dot{\gamma}^{-0.64}$",size=15)
        #plt.text(16000,20,"$C_3$="+str(2930))
    
        
        'SWITCH'
        #plt.xlim([0.1, 100000])
        #plt.xlim([0.1, 10])
        plt.xlim([100, 100000])
        plt.ylim([0.1, 100])
    
        'LABELS'
        #plt.xlabel("Units",size=13)
        #plt.xlabel("Alginate Volume Fraction, $\phi$ ($\%$)",size=13)
        plt.xlabel(r'Shear Rate, $\dot{\gamma}$ ($1/s$)',size=13)
        plt.ylabel(r'Microsphere Size, $d$ ($\mathrm{\mu m}$)',size=13)
        'log scale'
        plt.xscale('log')
        plt.yscale('log')
        plt.show()       
        
    if study=='f4':
        
        """TO CALCULATE  DIAMETER (FROM TAYLOR MODEL)""" 
        #alginate volume fraction
        phi = 100*np.array([0.00494364247577615,0.00739717920899497,0.00983864620228257,0.0122681322995387,0.0146857254748385,0.0170915128430511,0.019485580670304])
        
        #Alginate in dispersed phase, Algd
        #Calg=C/2
        
        
        'LENG AND CALABRESE MODEL'
        
        '''DISP. PHASE VISCOSITY'''
        phi_t=np.linspace(0.4,2)
        
        slengd1=lengTKE(phi=phi_t,S=9387)
        
        """DIAMETER (EXPERIMENTAL DATA) MEANS AND STDEVs"""         
        #mean (mu_d)
        d1=meanVmudx
        e1=sigmaVmudx
        
        '''MEDIANS'''            
        medmud=np.zeros(7)
        for i in range(0,7):
            medmud[i]=np.median(mudxM[i])
        
        boxwidth=0.001*np.ones(7)
        
        'ALGINATE FRACTION PLOT'
        
        CIl=np.array([0.104865625192891,0.292675903381101,0.464336520848279,0.624731530990819,0.776553584104539,0.921494853460006,1.06071046611521])
        CIu=np.array([0.961606244324977,1.14941652251319,1.32107713998037,1.48147215012291,1.63329420323662,1.77823547259209,1.9174510852473])
        #CIl=CIl[::-1]
        
        #plt.plot(phi_t,slengd1,color='magenta',label=r'Alginate Load, $\phi$ ($\%$)')
        plt.plot(phi_t,slengd1,color='k')
        plt.plot(phi,medmud,marker='o',color='r',linestyle='None')   
        'confidence intervals'
        plt.plot(phi,CIl,color='magenta')
        plt.plot(phi,CIu,color='magenta')
                
        plt.boxplot(mudxM,positions=phi,widths = 0.001)
        plt.text(0.13,3.5,"$\dot{\gamma}$="+str(9390)+r" $1/s$")
        #plt.text(0.13,6,"$d = e^{10.97}\phi^{0.92}\dot{\gamma}^{-0.74}$",size=15)
        plt.text(0.13,6,"$d = e^{9.22}\phi^{0.75}\dot{\gamma}^{-0.64}$",size=15)
        
        
        'SWITCH'
        #plt.xlim([0.1, 100000])
        plt.xlim([0.1, 10])
        'LABELS'
        #plt.xlabel("Units",size=13)
        plt.xlabel("Alginate Volume Fraction, $\phi$ ($\%$)",size=13)
        #plt.xlabel(r'Shear Rate, $\dot{\gamma}$ ($1/s$)',size=13)
        plt.ylabel(r'Microsphere Size, $d$ ($\mathrm{\mu m}$)',size=13)
        'log scale'
        plt.xscale('log')
        plt.yscale('log')
        plt.show()     
            
        
    if study=='f5':
        
        '''SHEAR RATES ARRAY'''
        RSc=3.82913349489977 
        O8K_gam= shear(N=8000,Ri=4,R=RSc)
        O6K_gam= shear(N=6000,Ri=4,R=RSc)
        O4K_gam= shear(N=4000,Ri=4,R=RSc)
        O2K_gam= shear(N=2000,Ri=4,R=RSc)
        S4_gam= shear(N=21200,Ri=2,R=1.77)
        S3_gam= shear(N=16800,Ri=2,R=1.77)
        S2_gam= shear(N=12400,Ri=2,R=1.77)
        S1_gam= shear(N=8000,Ri=2,R=1.77)
        gammaV = np.array([O8K_gam,S4_gam,O6K_gam,S3_gam,S2_gam,O4K_gam,S1_gam,O2K_gam])
        
        '''HISTOGRAMS'''
        bins =20
        # A= SUM(a_ij)
        histM = [[0 for j in range(8)] for i in range(np.shape(micsMoct[2])[0])]
        binsM = [[0 for j in range(8)] for i in range(np.shape(micsMoct[2])[0])]
        # Create histograms and normalize total count to 1 (data):
        for i in range(0,8):
            histM[i],binsM[i]= np.histogram(micsMoct[i], bins = bins)
            histM[i] = [ float(n)/sum(histM[i]) for n in histM[i]]
        # Set histograms' parameters (data)
        center = [[0 for j in range(8)] for i in range(20)]
        width=np.zeros(8)
        for i in range(0,8):
            center[i] = (binsM[i][:-1]+binsM[i][1:])/2
            width[i]= 1*(binsM[i][1]-binsM[i][0]) 
                
        f, ax = plt.subplots(8)            
        
        'median'
        medv=np.zeros(8)
        for i in range(0,8):
            medv[i]=np.median(micsMoct[i])
        
        'sample size'
        Nsize=np.zeros(8)
        for j in range(8):
            Nsize[j]=np.size(micsMoct[j])
        for k in range(0,8):
            ax[k].text(0.11,0.22,"$N$="+str(Nsize[k]))
        
        
        # Generate histograms + hist. annotations/description + hist. color style         
        cmap = mpl.cm.cool
        for k in range(0,8):
            #ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color=cmap((k+2) / float(11)))
            ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color='w')
            #ax[k].text(10,0.22,"$d_{32}$="+str(round(medv[k], 1))+r"$\; \mathrm{\mu m}$")
            ax[k].text(10,0.12,"$d$="+str(round(medv[k], 1)))
            ax[k].set_ylabel(''), ax[k].set_xscale('log'), ax[k].set_xlim([0.1, 100]),ax[k].set_ylim([0, 0.3])
        
        #diameter sig figs
        #for k in range(0,8):
        #    ax[k].text(10, 0.12,r"$d$="+str(sigfig.round_sig(meanVoct[k], 1))+r"$\pm$"+str(sigfig.round_sig(sigmaVoct[k], 1))+r"$\; \mathrm{\mu m}$")
        #Shear rate sig figs
        for k in range(0,4):
            ax[k].text(10,0.22,"$\dot{\gamma}$="+str(sigfig.round_sig(gammaV[k],3))+r" $1/s$")
        for k in range(4,8):
            ax[k].text(10,0.22,"$\dot{\gamma}$="+str(sigfig.round_sig(gammaV[k],2))+r" $1/s$")
        ax[7].set_xlabel('Microsphere diameter, $\; \mathrm{\mu m}$')
        
        # Fine-tune figure
        f.subplots_adjust(hspace=0,wspace=0,left=0.08,right=0.95)
        plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
        plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
        plt.setp(ax[7].get_xticklabels(), visible=True)
        plt.suptitle('Shear Rate, $\dot{\gamma}$ ($1/s$) \nContinuous Phase: 1-Octadecene',size=13)
        #plt.tsubplots_adjust(bottom=0.2)
        plt.show()
        
        print gammaV
        print d32oct
        
        
    if study=='f6':
    
        iron16=([1.127038,0.604433,0.811078,0.813274,0.985175,0.9025,0.946864,0.570791,0.815464,0.627596,0.655399,0.54522,0.517529,0.695022,0.759404,0.904803,0.578293,0.563717,0.760382,0.751139,0.746375,1.051724,0.714426,0.519821,0.682724,0.936764,0.902665,0.971655,0.864823,0.589746,0.393164,0.435153,0.578036,0.781778,0.742181,0.798706,1.166822,0.592512,0.502668,0.592512,0.648102,1.047334,0.969358,0.860861,0.755676,0.896054,0.784246,0.768935,1.073544,0.81783,0.865682,0.5,0.830455,0.621647,0.903323,0.579321,0.807957,0.967055,0.707737,0.807957,0.91963,0.759404,0.672414,0.735138,0.993588,1.072852,0.665973,0.748562,0.948276,0.612251,0.796096,0.627596,0.797029,1.087301,0.493116,0.810528,0.755676,0.717125,0.736552,0.976081,0.74158,0.483066,0.903488,0.482759,0.781208,0.802049,0.756266,0.469016,0.851836,0.928477,0.953746,0.823445,0.589746,0.57313,0.817467,0.490698,0.624747,0.953746,0.541116])
        iron50=([0.900585137561296,0.874766806861024,0.717209018399104,1.11817763925502,0.579771713530492,0.593874863832138,0.469329778768813,0.380984655989987,0.545837020976068,0.642284681814998,0.580868728160518,0.674200650749058,0.918087456953137,0.817588381146626,0.390882009522336,0.355036216362192,0.968041478927985,0.32703535245865,0.73301233814871,0.883456756679496,0.58849057440877,0.503363357230438,0.578672619246278,0.802657591776215,0.907626567175978,0.905519887583616,0.472034871941315,0.819921020355111,0.680778373836595,0.725153729891547,0.521989715007226,0.472034871941315,0.372536052451481,0.360375406434716,0.372536052451481,0.759458697866466,0.473381621816256,0.472034871941315,0.640299253009655,0.330905728036285,0.549324832956088,0.472034871941315,0.525635787601577,0.720750812869432,0.411510460923871,0.576468144111951,0.465242649168128,0.601331880665566,0.79227968265912,0.465242649168128,0.53165701249133,0.586323014283504,0.387610972856483,0.466609003502625,0.810550186653081,0.66373253836042,0.568685189153658,0.481383004625007,0.629267429963315,0.292073705590311,0.499554752522776,0.851159973773619,0.462497830822489,0.569803548521301,0.551637789851007,0.662772693261899,0.749332128784164,0.490553028687976,0.585236223702668,0.585236223702668,0.577571433435391,0.420690262209844,0.824566502405098,0.370823233942262,0.707377650962284,0.883456756679496,0.487950606993654,0.813685789025644,0.746779030633561,0.809764388904951,0.945417452982588,0.865989660359182,0.945417452982588,0.602389633252035,0.537610804071944,0.829185958731205,0.518317994998359,0.724275287740019,0.654070704917301,0.602389633252035,1.0464157913909,0.746779030633561,0.561928293552051,0.589571366089551,0.563060073739078,0.604499585888189,0.520768684761852,0.358604509199552,0.358604509199552,1.04458905309123,1.00355990348634,0.646237240240094])
        iron100=([0.709175309899033,0.749332128784164,0.94136853736232,0.840623546736597,0.606602199491942,0.972634043069759,0.866724484131922,0.868925226948759,0.749332128784164,0.766965885595732,0.783391001073122,0.648204481442857,0.766965885595732,0.59172702727032,0.625207597709509,0.542326778643488,0.953463704050084,1.12725022317394,0.601331880665566,0.548164694738268,1.20794557752678,1.24888688130694,0.86819226585084,0.820697093430275,0.990791445759808,0.682646082075529,0.473381621816256,0.81914421201363,1.29000933522547,0.767795486579845,0.942044565124394,0.852654549066489,0.925683646613013,0.942044565124394,0.830720069124647,1.07403679058497,0.927058084855655,0.851159973773619,1.03109408874935,1.02924015721291,0.907626567175978,0.95479815355413,1.41205280604348,1.02924015721291,0.884896784948462,1.23143343213492,1.3355929777346,0.832251351657616,1.00735888200478,0.925683646613013,1.13007046832017,0.802657591776215,0.90270333367641,0.802657591776215,1.12894321570214,0.648204481442857,0.641292735768508,0.953463704050084,0.984991641126052,0.726907429501923,0.625207597709509,0.781764019044672,0.711863286596938,0.648204481442857,0.888486644658096,0.925683646613013,1.01491417947073,0.725153729891547,1.30473034428993,1.13119659761636,0.805033495123152,0.884896784948462,1.05791188660597,0.890633615113424,0.682646082075529,1.0451983205989,0.628254931431422,0.833779821910676,1.20741843563934,0.847412003237287,0.728656908396924,1.0464157913909,0.783391001073122,0.889918533813782,1.02924015721291,0.942720108104809,1.04458905309123,0.950119380508171,0.889918533813782,0.834543007262139,0.767795486579845,0.605551805284838,0.576468144111951,0.925683646613013,0.905519887583616,0.851159973773619,0.490553028687976,0.817588381146626,0.884896784948462,0.628254931431422])
        iron200=([0.79067099128839,0.741646467884524,0.715431545980142,0.700140860629515,0.66946270975072,0.628254931431422,0.858606834588253,0.920856957835265,1.03294469285053,0.833015937350939,0.64525137049711,1.61362166228128,0.828417838178448,0.999746489174682,0.935942736911715,0.877673016883152,0.653096660140197,0.608697550116585,0.510895396995029,0.994639194009175,0.89917023466462,0.976553336488405,0.692828484089068,0.756939756606048,1.02985850521535,0.768624192149268,0.777681672504375,0.529256742840123,0.751029372473326,0.785014611107219,0.709175309899033,0.785014611107219,0.639304226374257,0.8368283871884,1.58736922580441,0.832251351657616,1.02676304132975,1.39070157418934,1.04397943001445,1.39070157418934,0.911126875450025,0.657952464247954,1.17752191666083,0.679842595560814,1.08816946134492,0.990791445759808,0.89917023466462,1.53393177119146,1.1710162209294,1.04397943001445,0.8368283871884,1.13063367317363,0.612866760150094,0.832251351657616,0.990791445759808,1.09458595441189,0.911126875450025,0.7882517984951,0.710968423532198,0.990791445759808,1.12725022317394,1.31057244967523,0.563060073739078,0.760296490396848,1.09050709983138,1.12328998430171,1.30619332028039,1.27811045527967,1.27011594084141,1.02303612144064,1.1672044863119,1.22105016558898,1.12328998430171,1.46645886101782,0.854891512811996,1.09050709983138,1.1672044863119,1.35170216641646,1.12328998430171,1.73491737087215,1.00165501109479,1.1672044863119,1.22677183774203,1.23864993093524,1.51429997690947,1.31154360390237,0.950119380508171,0.987573537144012,0.942720108104809,1.09050709983138,1.53807642011635,1.30619332028039,1.58174451331579,1.02303612144064,1.02303612144064,1.12328998430171,1.12328998430171,0.558519192562006,0.596014960366861,0.572033707920204])
        iron300=([0.639304226374257,0.751029372473326,0.639304226374257,0.735613218010963,0.639304226374257,0.762804347092994,0.958126146054919,0.783391001073122,0.73214332499416,0.750181230618936,0.882014377339272,0.73214332499416,1.02303612144064,0.67514424904252,0.762804347092994,0.882014377339272,0.510895396995029,1.0488464933684,0.934581364937436,0.826879456464352,0.880569635374576,1.23916378697515,0.717209018399104,0.880569635374576,0.840623546736597,0.952127384256663,1.02614282801956,0.60765077797465,0.952127384256663,0.974595659939255,0.618038723237103,0.699230995578931,0.916004836221115,0.880569635374576,0.783391001073122,0.952127384256663,0.751029372473326,0.597082132144185,1.30668061495143,0.39250730555361,0.767795486579845,0.840623546736597,0.840623546736597,0.976553336488405,0.766965885595732,0.844401649128951,0.829185958731205,0.738204934382485,0.829185958731205,0.785014611107219,0.999746489174682,0.677027500257308,0.828417838178448,0.825338207302505,0.724275287740019,0.764471680988537,0.73301233814871,0.867458685436074,0.6221453380035,0.666603781201826,0.828417838178448,0.780133643924814,0.6221453380035,0.954131162097561,0.718096104722579,0.768624192149268,0.931169227350564,0.751029372473326,0.666603781201826,0.867458685436074,0.768624192149268,0.780133643924814,0.692828484089068,0.735613218010963,0.807402407046451,1.01679423204744,1.29640896288571,0.735613218010963,0.798682046819832,0.6221453380035,0.519544784868306,0.545837020976068,0.766965885595732,1.14183943433777,1.02800234539138,0.867458685436074,0.867458685436074,0.677027500257308,0.677027500257308,0.735613218010963,0.785014611107219,0.653096660140197,0.807402407046451,0.735613218010963,0.751029372473326,0.931169227350564,0.702863388829752,0.808190501335631,0.865989660359182,0.717209018399104])
        
        MNP=[iron16,iron50,iron100,iron200,iron300]    
        
        '''HISTOGRAMS'''
        bins =20
        # A= SUM(a_ij)
        histM = [[0 for j in range(5)] for i in range(np.shape(MNP[2])[0])]
        binsM = [[0 for j in range(5)] for i in range(np.shape(MNP[2])[0])]
        # Create histograms and normalize total count to 1 (data):
        for i in range(0,5):
            histM[i],binsM[i]= np.histogram(MNP[i], bins = bins)
            histM[i] = [ float(n)/sum(histM[i]) for n in histM[i]]
        # Set histograms' parameters (data)
        center = [[0 for j in range(5)] for i in range(20)]
        width=np.zeros(5)
        for i in range(0,5):
            center[i] = (binsM[i][:-1]+binsM[i][1:])/2
            width[i]= 1*(binsM[i][1]-binsM[i][0]) 
            
        'median'
        medv=np.zeros(5)
        for i in range(0,5):
            medv[i]=np.median(MNP[i])  
                
        Zmnp=[16,50,100,200,300]
        
        f, ax = plt.subplots(5)    
    
        'sample size'
        Nsize=np.zeros(5)
        for j in range(5):
            Nsize[j]=np.size(MNP[j])
        for k in range(0,5):
            ax[k].text(0.11,0.22,"$N$="+str(Nsize[k]))
    
        for k in range(0,5):
            #ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color=cmap((2*k) / float(11)))
            ax[k].bar(center[k], histM[k], align = 'center', width = width[k],color='w')
            ax[k].text(10, 0.12,r"$d$="+str(round(medv[k], 1)))
            ax[k].set_ylabel(''), ax[k].set_xscale('log'), ax[k].set_xlim([0.1, 100]),ax[k].set_ylim([0, 0.3])
        
        for k in range(0,5):
            ax[k].text(10,0.22,"Z="+str(sigfig.round_sig(Zmnp[k],3))+r" mg/mL")
            ax[k].text(0.11,0.22,"$N$="+str(Nsize[k]))
            #ax[k].text(0.11,0.20,"$\mu_d$="+str(sigfig.round_sig(mudV[k],3))+r"$\; \mathrm{mPa-s}$")
            #ax[k].text(0.11,0.20,"$C_{alg,d}$="+str(sigfig.round_sig(Calg[k],2))+r"$\; \mathrm{mg/mL}$")
            
        ax[4].set_xlabel('Microsphere Size, $d$ ($\mathrm{\mu m}$)',size=13)
        
        # Fine-tune figure
        f.subplots_adjust(hspace=0,wspace=0,left=0.08,right=0.95)
        plt.setp([a.get_xticklabels() for a in f.axes[:]], visible=False)
        plt.setp([a.get_yticklabels() for a in f.axes[:]], visible=False)
        plt.setp(ax[4].get_xticklabels(), visible=True)
        #plt.suptitle(r"Alginate in Dispersed Phase, $C_{alg,d}$ ($\mathrm{mg/mL}$)",size=13)
        plt.suptitle(r"Magnetic Nanoparticle Concentration, Z (mg/mL)",size=13)
        #plt.tsubplots_adjust(bottom=0.2)
        #plt.show() 