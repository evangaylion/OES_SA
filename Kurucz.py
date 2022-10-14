# -*- coding: utf-8 -*-
'''
This is the parsing code for using the Kurucz database
Unlink NISTASD.py, the output from the Kurucz database is to be copied into a .txt file
and the Kurucz reader calls from that file
'''
import numpy as np
def kurucz_read(ion):
    wavl=[]
    Akil=[]
    code=[]
    Eil=[]
    Jil=[]
    Ekl=[]
    Jkl=[]
    f=open('{}.txt'.format(ion),'r')
    for line in f:
        row = line.split()
        wavl.append(row[0])
        Akil.append(row[1])
        Eil.append(row[4])
        Jil.append(row[5])
        Ekl.append(row[6])
        Jkl.append(row[7])
    wav=np.array([float(i) for i in wavl])
    Aki=np.array([float(i) for i in Akil])
    Ei=np.array([float(i) for i in Eil])
    Ek=np.array([float(i) for i in Ekl])
    Ji=np.array([float(i) for i in Jil])
    Jk=np.array([float(i) for i in Jkl])
    gi=2*Ji+1
    gk=2*Jk+1
    return wav, Aki, Ei, Ek, gi, gk
    
def kurucz(ion,lim_low,lim_hi):
    wav, Aki, Ei, Ek, gi, gk =kurucz_read(ion)
    i1 = np.argmin(np.abs(np.array(wav)-lim_low))
    i2 = np.argmin(np.abs(np.array(wav)-lim_hi))
    return wav[i1:i2],Aki[i1:i2],Ei[i1:i2],Ek[i1:i2],gi[i1:i2],gk[i1:i2]
