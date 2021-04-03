#!/usr/bin/env python3

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

# %%%


# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='1 for deterministic run, 0 for stochastic', default=1)
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)', default=96)
parser.add_argument('--Vn', metavar='Vn', help='the volume of the nucleus in liters', default=1.7500E-12)
parser.add_argument('--Vc', metavar='Vc', help='the volume of the cytoplasm in liters', default=5.2500E-12)
parser.add_argument('--folder', metavar='folder', help='input data folder path', default='input_files')
parser.add_argument('--outfile', metavar='outfile', help='the prefix for the name of the output files', default='output')
args = parser.parse_args()

input_data_folder = args.folder

if args.time == None or args.deterministic == None or args.Vn == None or args.Vc == None or args.outfile == None:
    print("ERROR: missing arguments. Need to pass --time, --deterministic, --Vn, --Vc, --outfile. Use -h for help.")

th = args.time
# th = 96
Vn = float(args.Vn)
Vc = float(args.Vc)
ts = 30
#
STIMligs = [100, 100.0, 100.0, 100.0, 100.0, 100.0, 1721.0]  # EGF, Her, HGF, PDGF, FGF, IGF, INS
# STIMligs = [100.0,0.0,0.0,0.0,0.0,0.0,100.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
# STIMligs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS

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



# %% Rate constants - kTL

# E2F 
model.setFixedParameterById('k9_1',0.63880304)
model.setFixedParameterById('k10_1',0.63880304 )
model.setFixedParameterById('k11_1',0.63880304 )
# model.setFixedParameterById('k12_1',0.1181584)
# model.setFixedParameterById('k13_1',0.1181584)
# model.setFixedParameterById('k14_1',0.1181584)


model.setFixedParameterById('k15_1',0.02*8) #CCNE1 
model.setFixedParameterById('k16_1',0.02*8) #CCNE2

# E2fatrep
model.setFixedParameterById('k150_1', model.getFixedParameterById('k150_1')*8)
model.setFixedParameterById('k151_1', model.getFixedParameterById('k151_1')*8)

# Ca
model.setFixedParameterById('k21_1', model.getFixedParameterById('k21_1')*100)

#%% Rate constants - CC

# model.setFixedParameterById('k328', 0.001) # pRB p
# model.setFixedParameterById('k329', 0.0001)
# model.setFixedParameterById('k330', 0.0001)
# model.setFixedParameterById('k331', 0.001) # pRB_E2F p

# model.setFixedParameterById('k332', 0.001)
# model.setFixedParameterById('k333', 0.0001)
# model.setFixedParameterById('k334', 0.001)
# model.setFixedParameterById('k335', 0.0001)
# model.setFixedParameterById('k336', 0.001)
# model.setFixedParameterById('k337', 0.001)
# model.setFixedParameterById('k338', 0.0001)
# model.setFixedParameterById('k339', 0.001)
# model.setFixedParameterById('k340', 0.001)

#%% Deterministic Run

flagD = 1
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints


xoutS_all, xoutG_all, xoutObs_all, tout_all = RunSPARCED(flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)



species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt', sep='\t', header=0, index_col=0).index
numberofgenes = int(len(genes_all))

obs_all = model.getObservableIds()

title = str()

for k in range(len(STIMligs)):
    title = title + STIMligs_id[k]+'='+str(STIMligs[k]) + " "


mpl.rcParams['figure.dpi'] = 300


#%% test module - diagnostics
from modules.RunSPARCED_test import RunSPARCED_test

flagD = 1
# nmxlsfile = outfile
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints


xoutS_all, xoutG_all, xoutObs_all, tout_all, vTC_all, vTCd_all, Nb_all, Nd_all, hills_all, TFa_all, TFr_all = RunSPARCED_test(flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)



species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt', sep='\t', header=0, index_col=0).index
numberofgenes = int(len(genes_all))

obs_all = model.getObservableIds()

title = str()

for k in range(len(STIMligs)):
    title = title + STIMligs_id[k]+'='+str(STIMligs[k]) + " "


mpl.rcParams['figure.dpi'] = 300


        
#%%
species_CC_dash_2 = ['ppERK', 'ppAKT', 'pcFos_cJun', 'cMyc', 'Cd', 'Cdk46',
                     'Cd_Cdk46', 'Cd_Cdk46_pRB', 'Cd_Cdk46_pRB_E2F', 'pRB',
                     'pRBp', 'pRBpp', 'pRB_E2F','pRBp_E2F','p107_E2Frep',
                     'p130_E2Frep', 'E2F', 'E2Frep', 'E2Fatrep', 'Cd_Cdk46_p18', 
                     'Cd_Cdk46_p19', 'Cd_Cdk46_p21', 'Cd_Cdk46_p27',
                     'Ce', 'Ce_Cdk2', 'Ce_Cdk2_pRBp', 'Ce_Cdk2_pRBp_E2F',
                     'Ce_Cdk2_p21', 'Ce_Cdk2_p57', 'Ce_Cdk2_p27',
                     'Ca', 'Ca_Cdk2', 'Ca_Cdk2_p21', 'Ca_Cdk2_p27',
                     'p18', 'p19', 'p21', 'p27', 'p107', 'p130']                     

k=0

cc_dash_species, axs_s = plt.subplots(10,4, sharex='col', figsize = (6,7))
plt.subplots_adjust(hspace = 0.6,wspace = 0.25)

cc_dash_species.suptitle(title,fontsize=5,y=0.92)

for i in range(10):
    for j in range(4):
        if k == len(species_CC_dash_2):
            break
        else:
            y_val = xoutS_all[:, list(species_all).index(species_CC_dash_2[k])]
            axs_s[i,j].plot(tout_all/3600, y_val, 'b-')            
            axs_s[i,j].set_ylim(0,max(y_val)*1.2)
            axs_s[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_s[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_s[i,j].title.set_text('species: '+species_CC_dash_2[k]+' (nM)')
            axs_s[i,j].title.set_size(5)
            if i == 9:
                axs_s[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1



#%%

obs_CC_dash = ['ERK', 'AKT', 'Fos', 'Jun', 'Myc', 'Cd', 'Cdk46', 'Cdk1', 'Cdk2', 'RB', 'E2F', 'Ce', 'Ca', 'Cb', 'Skp2', 'Pai', 'Pei', 'Pbi', 'p27', 'p107', 'p130', 'Cdc20','Cdh1a', 'Chk1', 'p21', 'p18', 'p19', 'p57', 'E2Frep','E2Fatrep']

k=0
obs_all = model.getObservableIds()
cc_dash_obs, axs_o = plt.subplots(10,3, sharex='col', figsize = (5,7))

cc_dash_obs.suptitle(title,fontsize=5,y=0.92)


plt.subplots_adjust(hspace = 0.8, wspace=0.25)

for i in range(10):
    for j in range(3):
        if k == len(obs_CC_dash):
            break
        else:
            y_val = xoutObs_all[:, list(obs_all).index(obs_CC_dash[k])]
            axs_o[i,j].plot(tout_all/3600, y_val, 'g-')            
            axs_o[i,j].set_ylim(0,max(y_val)*1.2)
            axs_o[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_o[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_o[i,j].title.set_text('obs: '+obs_CC_dash[k]+' (nM)')
            axs_o[i,j].title.set_size(5)
            if i == 9:
                axs_o[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1


#%% CDKi vs C

obs_Cd = xoutObs_all[:, list(obs_all).index('Cd')]
obs_Ce = xoutObs_all[:, list(obs_all).index('Ce')]
obs_Ca = xoutObs_all[:, list(obs_all).index('Ca')]
obs_Cb = xoutObs_all[:, list(obs_all).index('Cb')]

obs_C = obs_Cd+obs_Ce+obs_Ca+obs_Cb

obs_p18 = xoutObs_all[:, list(obs_all).index('p18')]
obs_p19 = xoutObs_all[:, list(obs_all).index('p19')]
obs_p21 = xoutObs_all[:, list(obs_all).index('p21')]
obs_p27 = xoutObs_all[:, list(obs_all).index('p27')]
obs_p57 = xoutObs_all[:, list(obs_all).index('p57')]

obs_Cdki = obs_p18 + obs_p19 + obs_p21 + obs_p27 + obs_p57

plt.plot(tout_all/3600, obs_C)
plt.plot(tout_all/3600, obs_Cdki)

#%%
#mrna_CC = list(genes_all[[5,6,7,8,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,141,142,143]])
# mrna_CC = list(genes_all[[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22,24,25,26,27,28,29,144,145,146,147,148]])
mrna_CC = list(genes_all[[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,25,26,27,28,29,144,145,146,147,148,149,150]])
#mrna_CC_data = pd.read_csv('input_files/OmicsData.txt', sep='\t', index_col=0, header=0)['Exp RNA'][[14,15,16,17,19,20,21,22,23,24,25,26,27,28,29]]


#x_m = xoutG_all[:,282:][:, list(genes_all).index(mrna_CC[2])]


cc_dash_mrna, axs_m = plt.subplots(9,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8)

cc_dash_mrna.suptitle(title,fontsize=5,y=0.92)


k=0
for i in range(9):
    for j in range(3):
        if k == len(mrna_CC):
            break
        else:
            y_val = xoutG_all[:, (numberofgenes*2):][:, list(genes_all).index(mrna_CC[k])]
            axs_m[i,j].plot(tout_all/3600, y_val,'r-')
            #axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
            axs_m[i,j].set_ylim(0,max(y_val)*1.2)
            axs_m[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_m[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_m[i,j].title.set_text('mrna: '+mrna_CC[k]+' (mpc)')
            axs_m[i,j].title.set_size(5)
            if i == 8:
                axs_m[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1


#%%
def timecourse(species, x_s=xoutS_all, tout_all=tout_all):

    x_t = x_s[:, list(species_all).index(species)]
    plt.scatter(tout_all/3600, x_t)
    plt.ylabel('species: '+str(species))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)

    plt.show


def timecourse_mrna(gene_symbol, x_g=xoutG_all[:, 282:], tout_all=tout_all):
    x_t = x_g[:, list(genes_all).index(gene_symbol)]
    plt.scatter(tout_all/3600, x_t)
    plt.ylabel('mrna_'+str(gene_symbol)+'_mpc')
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    plt.show

def timecourse_obs(obs_id, x_o = xoutObs_all, tout_all=tout_all):
    x_t = x_o[:, list(obs_all).index(obs_id)]
    plt.scatter(tout_all/3600, x_t)
    plt.ylabel('obs: '+str(obs_id))
    plt.xlabel('time(h)')
    plt.ylim(0, max(x_t)*1.25)
    plt.show



#%% stoichiometric matrix/observables check

stm = pd.read_csv('input_files/StoicMat.txt', sep='\t', header=0, index_col=0)

obsm = pd.read_csv('input_files/Observables.txt', sep='\t', header=0, index_col=0)



#%% test module with vTC output



#%% Dash vTC

mrna_CC = list(genes_all[[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,25,26,27,28,29,144,145,146,147,148,149,150]])



cc_dash_vTC, axs_vTC = plt.subplots(9,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8, wspace = 0.35)

cc_dash_vTC.suptitle(title,fontsize=5,y=0.92)


k=0
for i in range(9):
    for j in range(3):
        if k == len(mrna_CC):
            break
        else:
            y_val = vTC_all[:, list(genes_all).index(mrna_CC[k])]
            axs_vTC[i,j].plot(tout_all[1:]/3600, y_val,'r-')
            #axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
            axs_vTC[i,j].set_ylim(0,max(y_val)*1.2)
            axs_vTC[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_vTC[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_vTC[i,j].title.set_text('vTC: '+mrna_CC[k])
            axs_vTC[i,j].title.set_size(5)
            if i == 8:
                axs_vTC[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1

#%% Dash vTCd

cc_dash_vTCd, axs_vTCd = plt.subplots(9,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8, wspace = 0.35)

cc_dash_vTCd.suptitle(title,fontsize=5,y=0.92)


k=0
for i in range(9):
    for j in range(3):
        if k == len(mrna_CC):
            break
        else:
            y_val = vTCd_all[:, list(genes_all).index(mrna_CC[k])]
            axs_vTCd[i,j].plot(tout_all[1:]/3600, y_val,'r-')
            #axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
            axs_vTCd[i,j].set_ylim(0,max(y_val)*1.2)
            axs_vTCd[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_vTCd[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_vTCd[i,j].title.set_text('vTCd: '+mrna_CC[k])
            axs_vTCd[i,j].title.set_size(5)
            if i == 8:
                axs_vTC[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1

#%% Dash hills

cc_dash_hills, axs_hills = plt.subplots(9,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8, wspace = 0.35)

cc_dash_hills.suptitle(title,fontsize=5,y=0.92)


k=0
for i in range(9):
    for j in range(3):
        if k == len(mrna_CC):
            break
        else:
            y_val = hills_all[:, list(genes_all).index(mrna_CC[k])]
            axs_hills[i,j].plot(tout_all[1:]/3600, y_val,'r-')
            #axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
            axs_hills[i,j].set_ylim(0,max(y_val)*1.2)
            axs_hills[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_hills[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_hills[i,j].title.set_text('hills: '+mrna_CC[k])
            axs_hills[i,j].title.set_size(5)
            if i == 8:
                axs_hills[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1


#%% Dash vTCd vs repressors

vTCd_vs_rep, axs_vTCdrep = plt.subplots(3,1, sharex='col', figsize = (5,7))

e2f = ['E2F1', 'E2F2', 'E2F3']

for i in range(3):
    y_val = vTC_all[:,list(genes_all).index(e2f[i])]
    axs_vTCdrep[i].plot(xoutS_all[1:,list(species_all).index('E2Fatrep')], y_val, 'g-')
    axs_vTCdrep[i].set_ylim(0,max(y_val)*1.2)
    axs_vTCdrep[i].tick_params(axis='both', which='major', labelsize='4')
    axs_vTCdrep[i].ticklabel_format(useOffset=False, style='plain')
    axs_vTCdrep[i].title.set_text('vTCd '+str(e2f[i])+' vs E2Fatrep')
    axs_vTCdrep[i].title.set_size(5)
    axs_vTCdrep[i].set_xlabel('E2Fatrep')



