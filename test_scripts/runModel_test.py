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
parser.add_argument('--deterministic', metavar='flagD', type=int, help='0 for deterministic run, 1 for stochastic',default=1)
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)',default=12)
parser.add_argument('--Vn', metavar='Vn', help='the volume of the nucleus in liters',default=1.7500E-12)
parser.add_argument('--Vc', metavar='Vc', help='the volume of the cytoplasm in liters',default=5.2500E-12)
parser.add_argument('--folder', metavar='folder', help='input data folder path',default='input_files')
parser.add_argument('--outfile', metavar='outfile', help='the prefix for the name of the output files',default='output')
args = parser.parse_args()

input_data_folder = args.folder

if args.time == None or args.deterministic == None or args.Vn == None or args.Vc == None or args.outfile == None:
    print("ERROR: missing arguments. Need to pass --time, --deterministic, --Vn, --Vc, --outfile. Use -h for help.")

flagD = args.deterministic
th = args.time
Vn = float(args.Vn)
Vc = float(args.Vc)
outfile = args.outfile
ts = 30

STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS

species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_data_folder,'Species.txt'), encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)

species_initializations[155:162] = STIMligs

#%%

# reaction rate constant modification

slr = pd.read_csv('silenced.txt', header=None, sep='\t')
slr = [slr.values[i][0] for i in range(len(slr))]
slr = np.array(slr)

ratelaws = pd.read_csv('input_files/Ratelaws.txt',sep='\t',header=0,index_col=0)

slr_param = ['k'+str(list(ratelaws.index).index(slr[i])+1) for i in range(len(slr))]

#%%




sys.path.insert(0, os.path.abspath(model_output_dir))


model_module = importlib.import_module(model_name)
model = model_module.getModel()

model_param = np.array(model.getFixedParameterIds())

slr_param_actual = []
for i in range(len(model_param)):
    for j in range(len(slr_param)):
        if slr_param[j] == model_param[i] or slr_param[j]+'_1' == model_param[i]:
            slr_param_actual.append(model_param[i])
            
for i in range(len(slr_param_actual)):
    model.setFixedParameterById(slr_param_actual[i],0)


#%%

if flagD == 1:
    flagWr = 1
    nmxlsfile = outfile
    
    #sys.path.insert(0, os.path.abspath(model_output_dir))
    
    

    # species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_data_folder,'Species.txt'), encoding='latin-1')])

    # species_initializations = []
    # for row in species_sheet[1:]:
    #     species_initializations.append(float(row[2]))
    # species_initializations = np.array(species_initializations)

    #model_module = importlib.import_module(model_name)
    #model = model_module.getModel()
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    model.setTimepoints(np.linspace(0,ts)) # np.linspace(0, 30) # set timepoints
 

    xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model,input_data_folder)

    # if flagWr==1:
    #     columnsS=[ele for ele in model.getStateIds()]
    #     columnsG = columnsS[773:914]
    #     resa = [sub.replace('m_', 'ag_') for sub in columnsG]
    #     resi = [sub.replace('m_', 'ig_') for sub in columnsG]
    #     columnsG2 = np.concatenate((resa, resi), axis=None)
    #     condsSDF = pd.DataFrame(data=xoutS_all,columns=columnsS)
    #     condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
    #     condsSDF = None
    #     condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2)
    #     condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
    #     condsGDF = None

#%% stochastic test
flagD = 0


flagWr = 1
nmxlsfile = outfile

sys.path.insert(0, os.path.abspath(model_output_dir))
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_data_folder,'Species.txt'), encoding='latin-1')])

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

xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model,input_data_folder)

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

#%% test - plot trajectories
#xoutS_all[1960,:]

import matplotlib.pyplot as plt
species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt', sep='\t', header=0, index_col=0).index

import matplotlib as mpl
mpl.rcParams['figure.dpi']=300


def timecourse(species,x_s = xoutS_all, tout_all = tout_all):    
    
    x_t = x_s[:,list(species_all).index(species)]
    plt.scatter(tout_all/3600,x_t)
    plt.ylabel(str(species))
    plt.xlabel('time(h)')
    plt.ylim(0,max(x_t)*1.25)
    
    plt.show
 
def timecourse_mrna(gene_symbol,x_g = xoutG_all[:,282:], tout_all = tout_all):
    x_t = x_g[:,list(genes_all).index(gene_symbol)]
    plt.scatter(tout_all/3600,x_t)
    plt.ylabel('mrna_'+str(gene_symbol)+'_mpc')
    plt.xlabel('time(h)')
    plt.ylim(0,max(x_t)*1.25)
    plt.show
    
    
#%%

timecourse('pRB', x_s = xoutS_all[:1908,:],tout_all= tout_all[:1908])

#%%

timecourse('pRB', x_s = xoutS_all,tout_all= tout_all)

#xoutS_all[:,list(species_all).index('E')]

#%%

timecourse_mrna('CCND3')



#%%
elif flagD == 1:
    flagWr = 1
    nmxlsfile = outfile

    sys.path.insert(0, os.path.abspath(model_output_dir))
    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_data_folder,'Species.txt'), encoding='latin-1')])

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
