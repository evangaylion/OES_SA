# -*- coding: utf-8 -*-
import numpy as np
import NISTASD
import Kurucz
import GLV
import parsing_los_data2 as LOS
#various scientific constants 
heV=4.135667662e-15
keV=8.6173303e-5
eps=8.854187817e-12
q=1.60217662e-19
c=2.99792458e8
k=1.380658e-23
h=6.6260693e-34

#class that defines the wavelength, energy levels, statistical weights, and Einstein coefficient of each transition
#for a given range of wavelengths. Temperature and density are arbitrary variables here and are not used in later calculations
class Spectral:
    def __init__(self, ion, lim_low,lim_hi,temp,eden,condition):
        self.ion=ion
        self.lim_low=lim_low
        self.lim_hi=lim_hi
        self.temp=temp
        self.eden=eden
        self.condition=condition
    def NIST_spectral_dat(self):
        wav, Aki,Ei,Ek,gi,gk= NISTASD.process(self.ion, self.lim_low, self.lim_hi,self.temp,self.eden,self.condition)
        return np.array(wav), np.array(Aki),np.array(Ei),np.array(Ek),np.array(gi),np.array(gk)
    def Kurucz_spectral_dat(self):
        wav, Aki,Ei,Ek,gi,gk= Kurucz.kurucz(self.ion, self.lim_low, self.lim_hi)
        return np.array(wav), np.array(Aki),np.array(Ei),np.array(Ek),np.array(gi),np.array(gk)
        
#partition function
def partition_func(specdata_U,temp):
    wav, Aki,Ei,Ek,gi,gk=specdata_U[0],specdata_U[1],specdata_U[2],specdata_U[3],specdata_U[4], specdata_U[5]    
    U=sum(np.array(gk)*np.exp(-np.array(Ek)/(temp)))
    return U
    
#returns unbroadened spectral intensity for list of transitions
def emission_intensity(spec_data,ion, lim_low, lim_hi,temp,eden,condition):
    ION_P=Spectral(ion, 0,lim_hi,temp,eden,condition)
    Pspec_data=Spectral.Kurucz_spectral_dat(ION_P)
    Z=partition_func(Pspec_data,temp)
    wav, Aki,Ei,Ek,gi,gk=spec_data[0],spec_data[1],spec_data[2],spec_data[3],spec_data[4], spec_data[5]    
    wav_m=wav*1e-9
    N=eden*1e6
    #N=(gk/Z)*(eden*1e6)*np.exp(-Ek/temp)
    h=6.62607004e-34
    c=299792458

    I=(h*c*N/(4*np.pi*Z))*np.exp(-Ek/temp)*gk*Aki/wav_m
    return wav,I

#compute the values of FWHM for eah transition
def get_broadening(spec_data, eden, temp ,stark_width=.84e-10, instrum=.08, key='hwhm'):
    mass=9.10938356e-31
    wav,Ei,Ek=spec_data[0],spec_data[2],spec_data[3]
    stark = stark_width*(float(eden)/1e16) 
    difE=(Ek-Ei)
    natural_nu = difE/(4.135*1e-15)
    natural= natural_nu/299792458000000000
    doppler = np.array(wav)*np.sqrt((2*temp*1.602e-19)/mass)/c 
    FWHM_gauss   = doppler+instrum
    HWHM_gauss = FWHM_gauss/2.0
    FWHM_lorentz = stark+natural
    HWHM_lorentz = FWHM_lorentz/2.0
    if key=='hwhm':
        return [HWHM_gauss, HWHM_lorentz]
    if key=='fwhm':
        return [FWHM_gauss, FWHM_lorentz]

#apply broadening to emission lines
def emis_line_profiles(spec_data,ion,lim_low,lim_hi, eden, temp,condition=2):
    wav=spec_data[0]
    W=np.arange(lim_low,lim_hi,.1)
    I = emission_intensity(spec_data,ion, lim_low, lim_hi,temp,eden,condition)[1]
    HWHM_gauss = get_broadening(spec_data, eden, temp)[0]
    HWHM_lorentz=get_broadening(spec_data, eden, temp )[1]
    G_set, L_set, V_set=[],[],[]
    for i in range(len(wav)):
        X=W-wav[i]
        G=GLV.G(X, HWHM_gauss[i])*I[i]
        G_set.append(G)
        L=GLV.L(X , HWHM_lorentz[i])*I[i]
        L_set.append(L)
        V =GLV.V(X , HWHM_lorentz[i], HWHM_gauss[i])*I[i]
        V_set.append(V)
    return [W, G_set, L_set, V_set]        

#sum individual emission peaks over wavelength range, create continuous spectrum
def finalize_emission(spec_data,ion,lim_low,lim_hi,eden,temp,condition=2): 
    data=emis_line_profiles(spec_data,ion,lim_low,lim_hi, eden, temp)
    w,g,l,v=data[0],data[1],data[2],data[3]
    G=np.array(g).sum(axis=0)
    V=np.array(v).sum(axis=0)
    L=np.array(l).sum(axis=0)
    return [w, G, L, V]

#oscillator strength
def osc_strength(spec_data,eden,temp):
    wav, Aki,Ei,Ek,gi,gk=spec_data[0],spec_data[1],spec_data[2],spec_data[3],spec_data[4], spec_data[5]    
    w=wav*1e-9
    #v=c/(wav*10**-9)#Hz
    fo=1.5e4*Aki*w**2*gk/gi
    return np.array(fo)

#optical depth for computing self-absorption profile
def TAU(spec_data,eden,temp,L): #give l in nm
    import scipy
    ni=eden*1e-6
    wav, Aki,Ei,Ek,gi,gk=spec_data[0],spec_data[1],spec_data[2],spec_data[3],spec_data[4], spec_data[5]
    hg,hl=get_broadening(spec_data,eden,1 ,stark_width=.84e-10, instrum=.08, key='hwhm') #in nm
    VOIGTS=[]
    osc=osc_strength(spec_data,eden,temp)
    Ks=[]
    for i in range(len(wav)):
        W=np.arange(wav[i]-10,wav[i]+10,.1)
        X=W-wav[i]
        V=scipy.special.voigt_profile(X,hg[i],hl[i])
        Ks.append(eden*osc[i]*V*(wav[i]*1e-9)**2*8.817e-15)
    return np.array(Ks)*L

#compute unbroadened spectral intensity of self-absorption components        
def absorption_intensity(spec_data,ion, lim_low, lim_hi,temp,eden,L,condition):
    T=TAU(spec_data,eden,temp,L)
    inten=[]
    wav=spec_data[0]
    Io=emission_intensity(spec_data,ion, lim_low, lim_hi,temp,eden,condition)[-1]
    for i in range(len(wav)):
        integ=np.trapz((1-np.exp(-T[i])),axis=0)
        inten.append(Io[i]*integ)
    return wav,inten
   
#apply broadening to self-absorption component peaks
def abs_line_profiles(spec_data,ion, lim_low, lim_hi,temp,eden,L,condition ):
    wav=spec_data[0]
    W=np.arange(lim_low,lim_hi,.1)
    I = absorption_intensity(spec_data,ion, lim_low, lim_hi,temp,eden,L,condition)[1]
    HWHM_gauss = get_broadening(spec_data, eden, temp)[0]
    HWHM_lorentz=get_broadening(spec_data, eden, temp )[1]
    G_set, L_set, V_set=[],[],[]
    for i in range(len(wav)):
        X=W-wav[i]
        G=GLV.G(X, HWHM_gauss[i])*I[i]
        G_set.append(G)
        L=GLV.L(X , HWHM_lorentz[i])*I[i]
        L_set.append(L)
        V =GLV.V(X , HWHM_lorentz[i], HWHM_gauss[i])*I[i]
        V_set.append(V)
    return [W, G_set, L_set, V_set]        

#sum individual peaks over the full wavelength array to create continuous spectrum
def finalize_absorption(spec_data,ion, lim_low, lim_hi,eden,temp,L,condition): 
    data=abs_line_profiles(spec_data,ion, lim_low, lim_hi,temp,eden,L,condition)
    w,g,l,v=data[0],data[1],data[2],data[3]
    W=np.arange(lim_low,lim_hi,.1)
    G=np.array(g).sum(axis=0)
    V=np.array(v).sum(axis=0)
    L=np.array(l).sum(axis=0)
    return [w, G, L, V]

Fields='U00' #global variable for LOS_LINES
#integrate signal over a line-of-sight
def LOS_lines(timestep,fields,spec_data,ion, lim_low,lim_hi,L, condition=2):
    Ne=np.array(LOS.clean_densities(fields)[timestep])*1e-6
    Te=np.array(LOS.clean_temperatures()[timestep])/11606.0
    Wav,Voigt_Em, Voigt_Ab=[],[],[]
    for i in range(len(Ne)):
        Wav.append(finalize_emission(spec_data,ion,lim_low,lim_hi,Ne[i],Te[i],condition=2)[0])
        Voigt_Em.append(finalize_emission(spec_data,ion,lim_low,lim_hi,Ne[i],Te[i],condition=2)[-1])
        Voigt_Ab.append(finalize_absorption(spec_data,ion,lim_low,lim_hi,Ne[i],Te[i],L,condition=2)[-1])
    w1,g1,l1,v1=finalize_emission(spec_data,ion,lim_low,lim_hi,Ne[0],Te[0],condition=2)
    wa1,ga1,la1,va1=finalize_absorption(spec_data,ion,lim_low,lim_hi,Ne[0],Te[0],1e-3,condition=2)
    wavelengths=Wav[0]
    VvvE=np.array(Voigt_Em).sum(axis=0)
    VvvA=np.array(Voigt_Ab).sum(axis=0)
    
    LINES=VvvE-VvvA
    return [wavelengths,LINES]

#degree of ionization from Boltzmann Relation
def ionization(timestep,fields,lim_low,lim_hi):
    Ne=np.array(LOS.clean_densities(fields)[timestep])*1e-6
    Te=np.array(LOS.clean_temperatures()[timestep])/11606.0
    U=Spectral('U+I',lim_low,lim_hi,1.0,1e20,2)
    U_data=Spectral.Kurucz_spectral_dat(U)
    UII=Spectral('U+II',lim_low,lim_hi,1.0,1e20,2)
    UII_data=Spectral.Kurucz_spectral_dat(U)
    ionization_ratio=[]
    for i in range(len(Ne)):
        Z_I=partition_func(U_data,Te[i])
        Z_II=partition_func(UII_data,Te[i])
        ionization_ratio.append(((2*np.pi*Te[i])**(3.0/2.0)/4.1357e-15 )*(2*Z_II/Z_I)*np.exp(-6.194/Te[i])/Ne[i])
    return ionization_ratio, Ne

#calculate the McWhirter Criterion
def mcwhirter(timestep,ion,lim_low,lim_hi):
    U=Spectral(ion,lim_low,lim_hi,1.0,1e20,2)
    Te=np.array(LOS.clean_temperatures()[timestep])/11606.0
    spec_data=Spectral.Kurucz_spectral_dat(U)
    MC=[]
    wav, Aki,Ei,Ek,gi,gk=spec_data[0],spec_data[1],spec_data[2],spec_data[3],spec_data[4], spec_data[5]
    for i in Te: 
        for j in range(len(Ek)):
            MC.append( 1.6e12*(i*1.160e4)**(1/2.0)*(Ek[j]-Ei[j])**3)
    return MC
