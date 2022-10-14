# -*- coding: utf-8 -*-
'''
The NIST.py is a parsing code that connects to the NIST ASD webpage relating to the
selected ion, wavelength range, temperature, and density. The temperature and density are 
required to establish the NIST ASD url but are arbitrary and not used in later calculations
'''
import requests
import numpy as np
def process(ion, lim_low,lim_hi,temp,eden,condition):
    handle='https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra={}&limits_type=0&\
low_w={}&upp_w={}&unit=1&submit=Retrieve+Data&de=0&lte_out=1&temp={}&eden={}&format=1&\
line_out=1&remove_js=on&en_unit=1&output=0&bibrefs=1&page_size=15&show_obs_wl=1&unc_out=1&order_out=0&\
max_low_enrg=&show_av={}&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&\
forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on'.format(ion,lim_low,lim_hi,temp,eden,condition)
    r=requests.get(handle,stream=True)
    nl=[]
    for line in r.iter_lines():
        if str.encode('|') in line:
            g=line.split(str.encode('|'))
            nl.append(g)
    data=np.array(nl[4:-1])
    ion,wav,Aki=[],[],[]
    Ei_k=[]
    gi_k=[]
        #get info
    for j in data:wav.append(j[0].replace(str.encode(' '),str.encode('')))
    for k in data:Aki.append(k[2].replace(str.encode(' '),str.encode('')))
    for l in data:Ei_k.append(l[4].replace(str.encode(' '),str.encode('')))
    for m in data:gi_k.append(m[11].replace(str.encode(' '),str.encode('')))
    while str.encode('') in wav:
        wav.remove(str.encode(''))   
    while str.encode('') in Aki:
        Aki.remove(str.encode(''))  
    while str.encode('') in Ei_k:
        Ei_k.remove(str.encode(''))  
    while str.encode('') in gi_k:
        gi_k.remove(str.encode(''))  
    Ei,Ekn=[],[]
    for i in Ei_k:
        i=i.split(str.encode('-'))
        Ei.append(float(i[0].replace(str.encode(' '),str.encode(''))))
        Ekn.append(i[-1].replace(str.encode(' '),str.encode('')))
    for n, i in enumerate(Ekn):
        if str.encode('[') in i:
            Ekn[n] = i[1:-2]
    Ek=[]      
    for i in Ekn:
        Ek.append(float(i))
    gi,gk=[],[]
    for i in gi_k:
        i=i.split(str.encode('-'))
        gi.append(float(i[0].replace(str.encode(' '),str.encode(''))))
        gk.append(float(i[-1].replace(str.encode(' '),str.encode(''))))
    W,A=[],[]
    for i in wav:
        W.append(float(i))
    for j in Aki:
        A.append(float(j))
    return np.array(W), np.array(A),np.array(Ei),np.array(Ek),np.array(gi),np.array(gk)
    
#a function for reading the partition function from the NIST ASD. This was not used 
#because the ASD does not contain any data for uranium 
def partition(ion,temp):
    link='https://physics.nist.gov/cgi-bin/ASD/energy1.pl?encodedlist=XXT2&de=0&spectrum={}&units=1&format=2&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&g_out=on&lande_out=on&perc_out=on&biblio=on&temp={}&submit=Retrieve+Data'.format(ion,temp)
    f = requests.get(link)
    text=f.text
    D=text[-7:-1]
    D=D.strip('=')
    D=D.strip(' ')
    Z=float(D)
    return Z
