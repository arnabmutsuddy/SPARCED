#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 20:02:47 2021

@author: arnab
"""
import re
import matplotlib as mpl

import argparse
import scipy.stats
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import libsbml
import importlib
import amici
import amici.plotting
import os
import sys

sys.path.append(os.getcwd()[0:os.getcwd().rfind('/')]+'/sparced/bin')

from modules.RunSPARCED import RunSPARCED
from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep


# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='1 for deterministic run, 0 for stochastic', default=1)
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)', default=1)
parser.add_argument('--Vn', metavar='Vn', help='the volume of the nucleus in liters', default=1.7500E-12)
parser.add_argument('--Vc', metavar='Vc', help='the volume of the cytoplasm in liters', default=5.2500E-12)
parser.add_argument('--folder', metavar='folder', help='input data folder path', default='input_files')
parser.add_argument('--outfile', metavar='outfile', help='the prefix for the name of the output files', default='output')
args = parser.parse_args()

input_data_folder = args.folder

if args.time == None or args.deterministic == None or args.Vn == None or args.Vc == None or args.outfile == None:
    print("ERROR: missing arguments. Need to pass --time, --deterministic, --Vn, --Vc, --outfile. Use -h for help.")

flagD = args.deterministic
th = args.time
# th = 96
Vn = float(args.Vn)
Vc = float(args.Vc)
outfile = args.outfile
ts = 30

STIMligs = [100, 100.0, 100.0, 100.0, 100.0, 100.0, 1721.0]  # EGF, Her, HGF, PDGF, FGF, IGF, INS


STIMligs_id = ['E', 'H', 'HGF', 'P', 'F', 'I', 'INS']

species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
    os.path.join(input_data_folder, 'Species.txt'), encoding='latin-1')])

species_all = [species_sheet[k][0] for k in range(1,len(species_sheet))]

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)



for k in range(len(STIMligs)):
    species_initializations[species_all.index(STIMligs_id[k])] = STIMligs[k]

# %% model import

sys.path.insert(0, os.path.abspath(model_output_dir))


model_module = importlib.import_module(model_name)
model = model_module.getModel()

model_param = np.array(model.getFixedParameterIds())

#%%

flagD = 1
nmxlsfile = outfile
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints


# xoutS_all, xoutG_all, xoutObs_all, tout_all = RunSPARCED(flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)

ts = 30 # time-step to update mRNA numbers
NSteps = int(th*3600/ts)
tout_all = np.arange(0,th*3600+1,30)    
mpc2nM_Vc = (1E9/(Vc*6.023E+23))

genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs, geneIDs = RunPrep(flagD,Vn,model,input_data_folder)

mrna_IDs_sge = ['mrna_'+x for x in geneIDs]

paramIds = np.array(model.getFixedParameterIds())
params = np.array(model.getFixedParameters())

param_p = re.compile('mrna_\w+')

mrnaIds = np.array(list(filter(param_p.match,paramIds)))

mrna_i = np.nonzero(np.isin(paramIds,mrnaIds))[0]

species_id = list(pd.read_csv(os.path.join(input_data_folder,'Species.txt'),header=0,index_col=0,sep="\t").index)


params[mrna_i] = np.dot(mExp_mpc,mpc2nM_Vc)
model.setFixedParameters(params)


spdata = species_initializations

xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
xoutS_all[0,:] = spdata # 24hr time point     


xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
xoutG_all[0,:] = genedata

xoutObs_all = np.zeros(shape=(NSteps+1,len(model.getObservableIds())))


solver = model.getSolver() # Create solver instance
solver.setMaxSteps = 1e10

model.setInitialStates(xoutS_all[0,:])
xoutObs_all[0,:] = amici.runAmiciSimulation(model, solver)['y'][0,:]


import time

t0 = time.time()

t1 = []
t2 = []




for qq in range(NSteps):
    t1_0 = time.time()
    genedata,xmN_nM,AllGenesVec = SGEmodule(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs,spIDs,mrna_IDs_sge)
    t1_i = time.time() - t1_0
    t1.append(t1_i)
    xmN_nM = xmN_nM.reindex(index=mrnaIds)
    params[mrna_i] = xmN_nM.values
    model.setFixedParameters(params)
    model.setInitialStates(xoutS_all[qq,:])
    t2_0 = time.time()
    rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
    t2_i = time.time() - t2_0
    t2.append(t2_i)
    xoutS_all[qq+1,:] = rdata['x'][-1,:]
    # xoutS_all[qq+1,:] = spdata
    
    xoutG_all[qq+1,:] = genedata
    xoutObs_all[qq+1,:] = rdata['y'][-1,:]
    if rdata['x'][-1,species_id.index('PARP')] < rdata['x'][-1,species_id.index('cPARP')]:
        print('Apoptosis happened')
        break
xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
tout_all = tout_all[0:len(xoutS_all)]

t_total = time.time() - t0

#%%

t1 = np.array(t1)
t2 = np.array(t2)

#%%

species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt', sep='\t', header=0, index_col=0).index
numberofgenes = int(len(genes_all))

obs_all = model.getObservableIds()

title = str()

for k in range(len(STIMligs)):
    title = title + STIMligs_id[k]+'='+str(STIMligs[k]) + " "


mpl.rcParams['figure.dpi'] = 300


species_CC_dash = ['ppERK', 'ppAKT', 'pcFos_cJun', 'cMyc', 'Cd', 'Cdk46', 'Cd_Cdk46', 'Cd_Cdk46_pRB', 'Cd_Cdk46_pRB_E2F', 'pRB', 'pRBp', 'pRBpp', 'pRB_E2F','pRBp_E2F','Ce_Cdk2_p21','Ce_Cdk2_p27', 'E2F', 'Ce', 'Ce_Cdk2', 'Ce_Cdk2_pRBp', 'Ce_Cdk2_pRBp_E2F', 'Cd_Cdk46_p18', 'Cd_Cdk46_p19', 'Cd_Cdk46_p21', 'Cd_Cdk46_p27', 'p18', 'p21', 'p27', 'p57', 'E2Frep']                     

k=0

cc_dash_species, axs_s = plt.subplots(10,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.6)

cc_dash_species.suptitle(title,fontsize=5,y=0.92)

for i in range(10):
    for j in range(3):
        if k == len(species_CC_dash):
            break
        else:
            y_val = xoutS_all[:, list(species_all).index(species_CC_dash[k])]
            axs_s[i,j].plot(tout_all/3600, y_val, 'b-')            
            axs_s[i,j].set_ylim(0,max(y_val)*1.2)
            axs_s[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_s[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_s[i,j].title.set_text('species: '+species_CC_dash[k]+' (nM)')
            axs_s[i,j].title.set_size(5)
            if i == 9:
                axs_s[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1