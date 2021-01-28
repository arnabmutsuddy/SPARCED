#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 22:49:06 2021

@author: arnab
"""
#!/usr/bin/env python3

import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats
import argparse

from modules.RunSPARCED import RunSPARCED

#%%

# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name



parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='0 for deterministic run, 1 for stochastic', default = 0)
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)', default = 96)
parser.add_argument('--Vn', metavar='Vn', help='the volume of the nucleus in liters', default = 1.7500E-12)
parser.add_argument('--Vc', metavar='Vc', help='the volume of the cytoplasm in liters', default = 5.2500E-12)
parser.add_argument('--folder', metavar='folder', help='input data folder path',default='input_files')
parser.add_argument('--outfile', metavar='outfile', help='the prefix for the name of the output files', default = 'output')
args = parser.parse_args()


if args.time == None or args.deterministic == None or args.Vn == None or args.Vc == None or args.outfile == None:
    print("ERROR: missing arguments. Need to pass --time, --deterministic, --Vn, --Vc, --outfile. Use -h for help.")

flagD = args.deterministic
th = args.time
Vn = float(args.Vn)
Vc = float(args.Vc)
outfile = args.outfile
ts = 30


if flagD == 0:
    flagWr = 1
    nmxlsfile = outfile
    
    sys.path.insert(0, os.path.abspath(model_output_dir))

    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))
    species_initializations = np.array(species_initializations)

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = columnsS[773:914]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
        condsGDF = None

elif flagD == 1:
    flagWr = 1
    nmxlsfile = outfile

    sys.path.insert(0, os.path.abspath(model_output_dir))
    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open('Species.txt', encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))

    species_initializations = np.array(species_initializations)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver()          # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)

    if flagWr==1:
        columnsS=[ele for ele in model.getStateIds()]
        columnsG = columnsS[773:914]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
        condsGDF = None

#%% plot trajectories

import matplotlib.pyplot as plt
species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index

import matplotlib as mpl
mpl.rcParams['figure.dpi']=300


def timecourse(species,x_s = xoutS_all, tout_all = tout_all):    
    
    x_t = x_s[:,list(species_all).index(species)]
    plt.scatter(tout_all/3600,x_t)
    plt.ylabel(str(species))
    plt.xlabel('time(h)')
    
    plt.show

#%%

#plt.scatter(range(np.shape(xoutS_all)[0]),xoutS_all[:,list(species_all).index('Ca')])

#plt.scatter(tout_all/3600,xoutS_all[:,list(species_all).index('Ca')])

timecourse('Ca')
timecourse('Cb')


#%%

timecourse('Ce')

#%%

# cemal's code

import libsbml
import importlib
import amici
import amici.plotting
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import scipy.stats

from gm import gm
from gm_Prep import gm_Prep


#%%

flagD = 1
th = 96
ts = 30
NSteps = th*3600/ts
NSteps = int(NSteps)

# SBML model we want to import
sbml_file = 'SPARCED_Brep.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name

sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()
solver = model.getSolver() # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints

kGsRead = pd.read_csv('OmicsData_Brep.txt',header=0,index_col=0,sep="\t")
gExp_mpc = np.float64(kGsRead.values[:,0])
mExp_mpc = np.float64(kGsRead.values[:,1])
kGin = np.float64(kGsRead.values[:,2])
kGac = np.float64(kGsRead.values[:,3])
kTCleak = np.float64(kGsRead.values[:,4])
kTCmaxs = np.float64(kGsRead.values[:,5])
kTCd = np.float64(kGsRead.values[:,6])
kmRNAs = kGsRead.values[:,7]

# Read-in the activators matrix and assign concentrations of activators
TAsRead = pd.read_csv('TAs.csv',header=0,index_col=0)
TAs0 = np.float64(TAsRead.values)
# Read-in the repressors matrix and assign concentrations of repressors
TRsRead = pd.read_csv('TRs.csv',header=0,index_col=0)
TRs0 = np.float64(TRsRead.values)

#%%

# turn off synthesis and degradation of cyclinD mRNA
numStocCells = 1
kTCleak[9:12] = 0.0 # turn off transcription
kTCmaxs[9:12] = 0.0 # turn off transcription
kTCd[9:12] = 0.0 # turn off mRNA degradation
m = [1.0, 10.0, 60.0]
namexlsfile = 'Fig3J_det_96hr_'    

startTime = datetime.now()
print(startTime)

Vn = 1.75E-12
Vc = 5.25E-12
spdata0 = pd.read_csv('Species_Brep.txt',header=0,index_col=0,sep="\t")
spdata = np.double(spdata0.values[:,1])
genedata0, GenePositionMatrix, AllGenesVec, xgac_mpc_D, xgin_mpc_D, xgac_mpc, xgin_mpc, kTCleak2 \
= gm_Prep(flagD, gExp_mpc, mExp_mpc, kGin, kGac, kTCleak, kTCmaxs, kTCd)
tout_all = np.arange(0,th*3600+1,30)

condsS = []
condsG = []
for nn in range(numStocCells): 
    for pp in range(len(m)):
        xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
        xoutS_all[0,:] = spdata       
        genedata = genedata0
        genedata[291:294] = genedata[291:294]*m[pp]
        xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
        xoutG_all[0,:] = genedata
        for qq in range(NSteps):
            genedata,AllGenesVec = gm(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kGin,kGac, \
                                       kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,TAs0,TRs0)
            for ii,kk in enumerate(kmRNAs):
                model.setFixedParameterById(kk,genedata[ii+282]*(1E9/(Vc*6.023E+23)))
            model.setInitialStates(xoutS_all[qq,:])
            rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
            xoutS_all[qq+1,:] = rdata['x'][-1,:]
            xoutG_all[qq+1,:] = genedata
            if rdata['x'][-1,103] < rdata['x'][-1,105]:
                print('Apoptosis happened')
                break
        xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
        condsS.append(xoutS_all)
        xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
        condsG.append(xoutG_all)    
        condsSDF = pd.DataFrame(data=xoutS_all,columns=[ele for ele in model.getStateIds()]) 
        condsSDF.to_csv(namexlsfile+'S_'+str(pp)+'.txt',sep="\t")    
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all) 
        condsGDF.to_csv(namexlsfile+'G_'+str(pp)+'.txt',sep="\t")
        condsGDF = None    
        print(datetime.now() - startTime)
print(datetime.now())

#%%

namexls1 = 'Fig3J_det_96hr_S_0'
namexls2 = 'Fig3J_det_96hr_S_1'
namexls3 = 'Fig3J_det_96hr_S_2'

spids = [45,49,58,68] # Cyclin D, E, A, B active-forms

condsS1 = []
condsS2 = []
condsS3 = []
condsDF = pd.read_csv(namexls1+'.txt',header=0,index_col=0,sep="\t")
condsDF = np.double(condsDF.values[:])
condsDF = condsDF[~np.all(condsDF == 0, axis=1)]
condsS1= condsDF
condsDF = None
condsDF = pd.read_csv(namexls2+'.txt',header=0,index_col=0,sep="\t")
condsDF = np.double(condsDF.values[:])
condsDF = condsDF[~np.all(condsDF == 0, axis=1)]
condsS2= condsDF
condsDF = None
condsDF = pd.read_csv(namexls3+'.txt',header=0,index_col=0,sep="\t")
condsDF = np.double(condsDF.values[:])
condsDF = condsDF[~np.all(condsDF == 0, axis=1)]
condsS3= condsDF
condsDF = None

tt1 = (np.arange(0,len(condsS1))*30.0/3600.0)
tt2 = (np.arange(0,len(condsS2))*30.0/3600.0)
tt3 = (np.arange(0,len(condsS3))*30.0/3600.0)
yy1 = np.array(condsS1[:,spids[0]])
yy2 = np.array(condsS2[:,105])
yy3 = np.array(condsS3[:,105])
    
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('Cyclins D, E, A, and B active levels')
ax1.plot(tt1, np.array(condsS1[:,spids[0]]),'b',linewidth=4)
ax1.plot(tt2, np.array(condsS2[:,spids[0]]),'r',linewidth=4)
ax1.plot(tt3, np.array(condsS3[:,spids[0]]),'tab:orange',linewidth=4)

ax2.plot(tt1, np.array(condsS1[:,spids[1]]),'b',linewidth=4)
ax2.plot(tt2, np.array(condsS2[:,spids[1]]),'r',linewidth=4)
ax2.plot(tt3, np.array(condsS3[:,spids[1]]),'tab:orange',linewidth=4)

ax3.plot(tt1, np.array(condsS1[:,spids[2]]),'b',linewidth=4)
ax3.plot(tt2, np.array(condsS2[:,spids[2]]),'r',linewidth=4)
ax3.plot(tt3, np.array(condsS3[:,spids[2]]),'tab:orange',linewidth=4)

ax4.plot(tt1, np.array(condsS1[:,spids[3]]),'b',linewidth=4)
ax4.plot(tt2, np.array(condsS2[:,spids[3]]),'r',linewidth=4)
ax4.plot(tt3, np.array(condsS3[:,spids[3]]),'tab:orange',linewidth=4)

plt.savefig('Fig3J_1.png')

#%%