# -*- coding: utf-8 -*-
import os
import numpy as np
#need to convert VisIt .curve format to .txt
def curve2txt(folder):
    for filename in os.listdir(folder):
        infilename = os.path.join(folder,filename)
        if not os.path.isfile(infilename): continue
        oldbase = os.path.splitext(filename)
        newname = infilename.replace('.curve', '.txt')
        output = os.rename(infilename, newname)

#each key corresponds to a list of files containing the LoS measurements at a given timestep
def make_file_list(key):
    k=[]
    for root, dirs, files in os.walk(os.getcwd()):
        for f in files:
            if f.find(key)>=0:
                k.append(f)
    return k

#removes the # curve line of text from top row
def cleanup_files(file_list):
    for file in file_list:
        with open(file, 'r') as fin:
            data = fin.read().splitlines(True)
        with open(file, 'w') as fout:
            fout.writelines(data[1:])

#collects data from each file containing either temperature, density, concentration, etc. 
#there are N files for each physical variable, where N is the number of timesteps
def data_per_timestep(filename):
    entries=[]
    x,y=[],[]
    if 'temp' in filename:
        key='temperature'
    elif 'dens' in filename:
        key='density'
    else: key='concentration'
    for line in open(filename,'r'):
        d=line.split(' ')
        x.append(float(d[0]))
        y.append(float(d[1][:-1]))
    #print (key)
    H=np.array(list(zip(x,y)))
    return H

#This was initially defined to take all molecular components, but now only needs to use the concentration of U
def multiply_massfraction(molecule):
    if molecule =='U00': mass=238.02891
    elif molecule =="UO00": mass=254.028
    elif molecule =="UO200":mass=270.027
    elif molecule =="UO300":mass=286.026
    elif molecule =="O00":mass=16.0
    elif molecule =="O200":mass=32.0
    else: print ("incorrect molecule entry")
    mass_components=[]
    #for molecule in fields:
    timesteps=make_file_list(molecule)
    
    for t in timesteps:
        mass_components.append(data_per_timestep(t)[0:-1][:,-1]*mass)
    return mass_components

#Convert the molecular concentration to number density
def collect_N(molecule='U00'):
    J=multiply_massfraction('U00')
    densities_f=make_file_list('density00')
    dens=[]
    for d in densities_f:
        dens.append(data_per_timestep(d)[:,-1]
    P=[] 
    I=[]
    for timestep in range(len(J)):
        k=[]
        denses=dens[timestep]
        indexx=[]
        for xdata in range(len(J[timestep])):
            XX=J[timestep][xdata]
            if XX>1e-4:
                number_density=denses[xdata]*1e3*XX*1e-1*6.0221409e+23
                k.append(number_density)
                indexx.append(xdata)
        P.append(k)
        I.append(indexx)
    return P,I

#returns only number density and temperature along LOS containing 100% uranium
def clean_densities(molecule):
    data,ind=collect_N(molecule)
    return data
def clean_temperatures(molecule='U00'):
    temperature_file_list= make_file_list('temp00') 
    indexes=collect_N(molecule='U00')[-1]
    D=[]
    for i in temperature_file_list:
        D.append(data_per_timestep(i)[:,-1])
    TT=[]   
    for Is in range(len(indexes)):
        temperature_list=D[Is]
        T_news=[]
        for Ksin in indexes[Is]:
            T_news.append(temperature_list[Ksin])
        TT.append(T_news)
    return TT
