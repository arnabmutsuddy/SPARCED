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

# %%


# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--deterministic', metavar='flagD', type=int, help='1 for deterministic run, 0 for stochastic', default=1)
parser.add_argument('--time', metavar='time', type=int, help='experiment run time (in hours)', default=48)
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
Vn = float(args.Vn)
Vc = float(args.Vc)
outfile = args.outfile
ts = 30

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

# %% turn off irrelevant CC reactions

#ks = ['k'+str(n) for n in range(320, 432)]

#ks = ['k'+str(n) for n in range(337, 432)] # turn on E2F, pRB reactions

# ks = ['k'+str(n) for n in range(328, 419)]

# ks_actual = []
# for p in model_param:
#     for k in ks:
#         if k in p:
#             ks_actual.append(p)

# for p in ks_actual:
#     model.setFixedParameterById(p, 0)

# %% Rate constants - kTL


# model.setFixedParameterById('k424', 0.001)  # kon, Mdi
# model.setFixedParameterById('k425', 0.0001)  # koff, Mdi
# model.setFixedParameterById('k426',0.001) #kon, Mei
# model.setFixedParameterById('k427',0.0001) #koff, Mei
# model.setFixedParameterById('k428',0.001) #kon, Mai
# model.setFixedParameterById('k429',0.0001) #koff, Mai
# model.setFixedParameterById('k430',0.001) #kon, Mbi
# model.setFixedParameterById('k431',0.0001) #koff, Mbi
# model.setFixedParameterById('k17_1', model.getFixedParameterById('k17_1')/12) #Skp2
# model.setFixedParameterById('k20_1', model.getFixedParameterById('k20_1')/25) #Pbi
# model.setFixedParameterById('k22_1', model.getFixedParameterById('k22_1')*1.5*1.25) #p27
# model.setFixedParameterById('k23_1', model.getFixedParameterById('k23_1')/15) #Cdh1a
# model.setFixedParameterById('k25_1', model.getFixedParameterById('k25_1')/6) #Cdc20
# model.setFixedParameterById('k26_1', model.getFixedParameterById('k26_1')/200) #Wee1
# model.setFixedParameterById('k28_1', model.getFixedParameterById('k28_1')*1.2*1.1) #p21


# model.setFixedParameterById('k31_1', model.getFixedParameterById('k31_1')*1.5) #Cdk46
# model.setFixedParameterById('k32_1', model.getFixedParameterById('k32_1')*1.5*1.1) #Cdk46



# model.setFixedParameterById('k12_1',model.getFixedParameterById('k12_1')*50)
# model.setFixedParameterById('k13_1',model.getFixedParameterById('k13_1')*50)
# model.setFixedParameterById('k14_1',model.getFixedParameterById('k14_1')*50)

# model.setFixedParameterById('k8_1',model.getFixedParameterById('k8_1')/450*1.1) #Rb
# model.setFixedParameterById('k9_1',model.getFixedParameterById('k9_1')/100/1.1)
# model.setFixedParameterById('k10_1',model.getFixedParameterById('k10_1')/100/1.1)
# model.setFixedParameterById('k11_1',model.getFixedParameterById('k11_1')/100/1.1)

# model.setFixedParameterById('k432', 1.47e-4)
# model.setFixedParameterById('k433', 1.47e-4)

# Cyclins

model.setFixedParameterById('k12_1',model.getFixedParameterById('k12_1')*5) #0.590792
model.setFixedParameterById('k13_1',model.getFixedParameterById('k13_1')*5)
model.setFixedParameterById('k14_1',model.getFixedParameterById('k14_1')*5)

model.setFixedParameterById('k15_1',model.getFixedParameterById('k15_1')*5*10*2) #CCNE1
model.setFixedParameterById('k16_1',model.getFixedParameterById('k16_1')*5*10*2) #CCNE2
# model.setFixedParameterById('k21_1',model.getFixedParameterById('k21_1')/100/50) #CCNA2
# model.setFixedParameterById('k24_1',model.getFixedParameterById('k24_1')/100/50) #CCNB2

#dParameterById('k144_1',0.14244776) #CDKN2C


#%% Rate constants - half lives

# # vTLCd
# model.setFixedParameterById('k154', 2.14e-5) # pRB half life
# model.setFixedParameterById('k155', 1.65e-4) # E2F half life


# # vCd
# model.setFixedParameterById('k341', 1.47e-4) # Cd_Cdk46, vCd1
# #model.setFixedParameterById('k433', 1.47e-4) # Md, vCd2 (obsolete)
# model.setFixedParameterById('k342', 1.47e-4) # Cd_Cdk46_p27, vCd3
# model.setFixedParameterById('k343', 3.85e-5) # Ce_Cdk2, vCd4
# #model.setFixedParameterById('k436', 3.85e-5) # Me, vCd5 (obsolete)
# model.setFixedParameterById('k344', 3.85e-5) # Ce_Cdk2_p27, vCd6
# model.setFixedParameterById('k345', 3.21e-5) # Pe, vCd7
# model.setFixedParameterById('k346', 7.70e-5) # Ca_Cdk2, vCd8
# # model.setFixedParameterById('k440', 7.70e-5) # Ma, vCd9 (obsolete)
# model.setFixedParameterById('k347', 7.70e-5) # Ca_Cdk2_p27, vCd10
# model.setFixedParameterById('k348', 2.35e-5) # Cdh1i, vCd11
# model.setFixedParameterById('k349', 1.65e-4) # E2Fp, vCd12
# model.setFixedParameterById('k350', 3.21e-5) # p27p, vCd13
# model.setFixedParameterById('k351', 3.85e-4) # Pa, vCd14
# model.setFixedParameterById('k352', 1.28e-4) # Cb_Cdk1, vCd15
# # model.setFixedParameterById('k353', 1.28e-4) # Mb, vCd16
# model.setFixedParameterById('k353', 3.85e-4) # Cdc20a, vCd17
# model.setFixedParameterById('k354', 1.93e-4) # Pb, vCd18
# model.setFixedParameterById('k355', 2.75e-5) # Wee1p, vCd19
# model.setFixedParameterById('k356', 1.28e-4) # Cb_Cdk1_p27, vCd20
# model.setFixedParameterById('k357', 1.65e-4) # pRB_E2F, vCd21
# model.setFixedParameterById('k358', 1.65e-4) # pRB_E2Fp, vCd22
# model.setFixedParameterById('k359', 3.85e-4) # Cd_Cdk46_p21, vCd23
# model.setFixedParameterById('k360', 3.85e-4) # Ce_Cdk2_p21, vCd24
# model.setFixedParameterById('k361', 3.85e-4) # Ca_Cdk2_p21, vCd25
# model.setFixedParameterById('k362', 3.85e-4) # Cb_Cdk1_p21, vCd26
# model.setFixedParameterById('k363', 1.65e-4) # pRBpp_E2F, vCd27
# model.setFixedParameterById('k364', 1.47e-4) # Cd_Cdk46_pRB, vCd28
# model.setFixedParameterById('k365', 1.65e-4) # Cd_Cdk46_pRB_E2F, vCd29
# model.setFixedParameterById('k366', 3.85e-5) # Ce_Cdk2_pRBp, vCd30
# model.setFixedParameterById('k367', 1.65e-4) # Ce_Cdk2_pRBp_E2F, vCd31

# # model.setFixedParameterById('k368', 1.47e-4) # Cd_Cdk46_pRBp, vCd32
# # model.setFixedParameterById('k369', 1.65e-4) # Cd_Cdk46_pRBp_E2F, vCd33
# # model.setFixedParameterById('k370', 3.85e-5) # Ce_Cdk2_pRBpp, vCd34
# # model.setFixedParameterById('k371', 1.65e-4) # Ce_Cdk2_pRBpp_E2F, vCd35

# model.setFixedParameterById('k368', 2.14e-5) # pRBp half life
# model.setFixedParameterById('k369', 2.14e-5) # pRBpp half life



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









#%% initial conditions
# species_all = [species_sheet[n][0] for n in range(1,len(species_sheet))]

# species_initializations[species_all.index('pRB')] = 21.484-0.2182*0.90
# species_initializations[species_all.index('pRB_E2F')] = 0.2182*0.90
# species_initializations[species_all.index('E2F')] = 0.2182*0.10

# species_initializations[species_all.index('p18')] = 40.7




#%% Deterministic Run
#model_module = importlib.import_module(model_name)
#model = model_module.getModel()
flagD = 1
nmxlsfile = outfile
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints


xoutS_all, xoutG_all, xoutObs_all, tout_all = RunSPARCED(flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)

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

#xoutS_final = pd.Series(data=xoutS_all[-1], index=model.getStateIds())

#xoutS_final.to_csv('xoutS_final.csv', sep='\t', index=True)
# %% test - plot trajectories
# xoutS_all[1960,:]

species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt', sep='\t', header=0, index_col=0).index
numberofgenes = int(len(genes_all))

obs_all = model.getObservableIds()

mpl.rcParams['figure.dpi'] = 300

 #%%#
species_CC_dash = ['ppERK', 'ppAKT', 'pcFos_cJun', 'cMyc', 'Cd', 'Cdk46', 'Cd_Cdk46', 'Cd_Cdk46_pRB', 'Cd_Cdk46_pRB_E2F', 'pRB', 'pRBp', 'pRBpp', 'pRB_E2F','pRBp_E2F','Ce_Cdk2_p21','Ce_Cdk2_p27', 'E2F', 'Ce', 'Ce_Cdk2', 'Ce_Cdk2_pRBp', 'Ce_Cdk2_pRBp_E2F', 'Cd_Cdk46_p18', 'Cd_Cdk46_p19', 'Cd_Cdk46_p21', 'Cd_Cdk46_p27', 'p18', 'p19', 'p21', 'p27', 'p57']                     

k=0

cc_dash_species, axs_s = plt.subplots(10,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.6)

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
        
#%%

obs_CC_dash = ['ERK', 'AKT', 'Fos', 'Jun', 'Myc', 'Cd', 'Cdk46', 'Cdk1', 'Cdk2', 'RB', 'E2F', 'Ce', 'Ca', 'Cb', 'Skp2', 'Pai', 'Pei', 'Pbi', 'p27', 'Cdh1a', 'Cdc20', 'Wee1', 'Chk1', 'p21', 'p18', 'p19', 'p57']

k=0
obs_all = model.getObservableIds()
cc_dash_obs, axs_o = plt.subplots(9,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8, wspace=0.25)

for i in range(9):
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
            if i == 8:
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
mrna_CC = list(genes_all[[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]])
#mrna_CC_data = pd.read_csv('input_files/OmicsData.txt', sep='\t', index_col=0, header=0)['Exp RNA'][[14,15,16,17,19,20,21,22,23,24,25,26,27,28,29]]


#x_m = xoutG_all[:,282:][:, list(genes_all).index(mrna_CC[2])]


cc_dash_mrna, axs_m = plt.subplots(9,3, sharex='col')
plt.subplots_adjust(hspace = 0.8)
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
# %%

timecourse('ppERK')
timecourse('ppAKT')

timecourse('pRB')
timecourse('E2F')


timecourse_obs('RB')
timecourse_obs('E2F')

timecourse_obs('Cd')
timecourse_obs('Cdk46')


timecourse('pRB')
timecourse('pRBp')
timecourse('pRBpp')
timecourse('pRB_E2F')
timecourse('pRBp_E2F')
timecourse('pRBpp_E2F')
timecourse('Cd_Cdk46_pRB')
timecourse('Cd_Cdk46_pRB_E2F')
timecourse('Ce_Cdk2_pRBp')
timecourse('Ce_Cdk2_pRBp_E2F')
timecourse('Cd_Cdk46_pRBp')
timecourse('Cd_Cdk46_pRBp_E2F')
timecourse('Ce_Cdk2_pRBpp')
timecourse('Cd_Cdk2_pRBpp_E2F')


timecourse('Cd')
timecourse('Cdk46')
timecourse('Cd_Cdk46')

timecourse('Cd_Cdk46_pRB')
timecourse('Cd_Cdk46_pRB_E2F')


timecourse('Cd_Cdk46_p27')
timecourse('Cd_Cdk46_p21')

timecourse('Ce_Cdk2')
timecourse('pRB')
timecourse('pRBp')
timecourse('pRBpp')
timecourse('E2F')
# xoutS_all[:,list(species_all).index('E')]
timecourse_obs('Cd')
timecourse_obs('RB')
timecourse_obs('E2F')

timecourse('p53ac')
# %%
timecourse_mrna('CCND1')
timecourse_mrna('CCND2')
timecourse_mrna('CCND3')

timecourse_mrna('RB1')
timecourse_mrna('E2F1')
timecourse_mrna('E2F2')
timecourse_mrna('E2F3')

timecourse_mrna('CCNE1')
timecourse_mrna('CCNE2')

timecourse_mrna('CCNA2')

#%% stoichiometric matrix/observables check

stm = pd.read_csv('input_files/StoicMat.txt', sep='\t', header=0, index_col=0)

obsm = pd.read_csv('input_files/Observables.txt', sep='\t', header=0, index_col=0)

#%%
species_CC_dash = ['ppERK', 'ppAKT', 'pcFos_cJun', 'cMyc', 'Cd', 'Cdk46', 'Cd_Cdk46', 'Cd_Cdk46_pRB', 'Cd_Cdk46_pRB_E2F', 'pRB', 'pRBp', 'pRBpp', 'pRB_E2F', 'E2F', 'Ce', 'Ce_Cdk2', 'Ce_Cdk2_pRBp', 'Ce_Cdk2_pRBp_E2F']

k=0

cc_dash_species, axs_s = plt.subplots(6,3, sharex='col')
plt.subplots_adjust(hspace = 0.6)

for i in range(6):
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
            if i == 5:
                axs_s[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1
        
#%%

obs_CC_dash = ['ERK', 'AKT', 'Fos', 'Jun', 'Myc', 'Cd', 'Cdk46', 'Cdk1', 'Cdk2', 'RB', 'E2F', 'Ce', 'Ca', 'Cb', 'Skp2', 'Pai', 'Pei', 'Pbi', 'p27', 'Cdh1a', 'Cdc20', 'Wee1', 'Chk1', 'p21',]

k=0
obs_all = model.getObservableIds()
cc_dash_obs, axs_o = plt.subplots(8,3, sharex='col', figsize = (5,7))
plt.subplots_adjust(hspace = 0.8, wspace=0.25)

for i in range(8):
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
            if i == 7:
                axs_o[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1


#%%
mrna_CC = list(genes_all[[5,6,7,8,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]])

#mrna_CC_data = pd.read_csv('input_files/OmicsData.txt', sep='\t', index_col=0, header=0)['Exp RNA'][[14,15,16,17,19,20,21,22,23,24,25,26,27,28,29]]


x_m = xoutG_all[:,282:][:, list(genes_all).index(mrna_CC[2])]


cc_dash_mrna, axs_m = plt.subplots(8,3, sharex='col')
plt.subplots_adjust(hspace = 0.8)
k=0
for i in range(8):
    for j in range(3):
        if k == len(mrna_CC):
            break
        else:
            y_val = xoutG_all[:, 282:][:, list(genes_all).index(mrna_CC[k])]
            axs_m[i,j].plot(tout_all/3600, y_val,'r-')
            #axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
            axs_m[i,j].set_ylim(0,max(y_val)*1.2)
            axs_m[i,j].tick_params(axis='both', which='major', labelsize='4')
            axs_m[i,j].ticklabel_format(useOffset=False, style='plain')
            axs_m[i,j].title.set_text('mrna: '+mrna_CC[k]+' (mpc)')
            axs_m[i,j].title.set_size(5)
            if i == 7:
                axs_m[i,j].set_xlabel('time(h)', fontsize=5)
            k +=1



#%% CC dashboard
species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt',sep='\t', header=0, index_col=0).index

obs_all = model.getObservableIds()


plt.yticks(fontsize=5)

cc_dash, axs = plt.subplots(4,4, sharex='col')

plt.subplots_adjust(hspace = 0.4)
axs[0,0].plot(tout_all/3600,xoutS_all[:, list(species_all).index('ppERK')], 'g-')
axs[0,0].tick_params(axis='both', which='major', labelsize='4')
axs[0,0].title.set_text('species: ppERK')
axs[0,0].title.set_size(5)
axs[0,1].plot(tout_all/3600,xoutS_all[:, list(species_all).index('ppAKT')], 'b-')
axs[0,1].tick_params(axis='both', which='major', labelsize='4')
axs[0,1].title.set_text('species: ppAKT')
axs[0,1].title.set_size(5)
axs[0,2].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pcFos_cJun')], 'c-')
axs[0,2].tick_params(axis='both', which='major', labelsize='4')
axs[0,2].title.set_text('species: pcFos_cJun')
axs[0,2].title.set_size(5)
axs[0,3].plot(tout_all/3600,xoutS_all[:, list(species_all).index('cMyc')], 'r-')
axs[0,3].tick_params(axis='both', which='major', labelsize='4')
axs[0,3].title.set_text('species: cMyc')
axs[0,3].title.set_size(5)
axs[1,0].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Cd')], 'r-')
axs[1,0].tick_params(axis='both', which='major', labelsize='4')
axs[1,0].title.set_text('species: Cd')
axs[1,0].title.set_size(5)
axs[1,1].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Cdk46')], 'm-')
axs[1,1].tick_params(axis='both', which='major', labelsize='4')
axs[1,1].title.set_text('species: Cdk46')
axs[1,1].title.set_size(5)
axs[1,2].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Cd_Cdk46')], 'y-')
axs[1,2].tick_params(axis='both', which='major', labelsize='4')
axs[1,2].title.set_text('species: Cd_Cdk46')
axs[1,2].title.set_size(5)
axs[1,3].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Cd_Cdk46_pRB')], 'k-')
axs[1,3].tick_params(axis='both', which='major', labelsize='4')
axs[1,3].title.set_text('species: Cd_Cdk46_pRB')
axs[1,3].title.set_size(5)
axs[2,0].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Cd_Cdk46_pRB_E2F')], 'g-')
axs[2,0].tick_params(axis='both', which='major', labelsize='4')
axs[2,0].title.set_text('species: Cd_Cdk46_pRB_E2F')
axs[2,0].title.set_size(5)
axs[2,1].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pRB')], 'b-')
axs[2,1].tick_params(axis='both', which='major', labelsize='4')
axs[2,1].title.set_text('species: pRB')
axs[2,1].title.set_size(5)
axs[2,2].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pRBp')], 'c-')
axs[2,2].tick_params(axis='both', which='major', labelsize='4')
axs[2,2].title.set_text('species: pRBp')
axs[2,2].title.set_size(5)
axs[2,3].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pRBpp')], 'r-')
axs[2,3].tick_params(axis='both', which='major', labelsize='4')
axs[2,3].title.set_text('species: pRBpp')
axs[2,3].title.set_size(5)
axs[3,0].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pRB_E2F')], 'm-')
axs[3,0].tick_params(axis='both', which='major', labelsize='4')
axs[3,0].title.set_text('species: pRB_E2F')
axs[3,0].title.set_size(5)
axs[3,0].set_xlabel('time(h)',fontsize=5)
axs[3,1].plot(tout_all/3600,xoutS_all[:, list(species_all).index('pRBp_E2F')], 'y-')
axs[3,1].tick_params(axis='both', which='major', labelsize='4')
axs[3,1].title.set_text('species: pRBp_E2F')
axs[3,1].title.set_size(5)
axs[3,2].plot(tout_all/3600,xoutS_all[:, list(species_all).index('E2F')], 'k-')
axs[3,2].tick_params(axis='both', which='major', labelsize='4')
axs[3,2].title.set_text('species: E2F')
axs[3,2].title.set_size(5)
axs[3,3].plot(tout_all/3600,xoutS_all[:, list(species_all).index('Ce')], 'g-')
axs[3,3].tick_params(axis='both', which='major', labelsize='4')
axs[3,3].title.set_text('species: Ce')
axs[3,3].title.set_size(5)

#%% plot - mrnas



# def timecourse_mrna(gene_symbol, x_g=xoutG_all[:, 282:], tout_all=tout_all):
#     x_t = x_g[:, list(genes_all).index(gene_symbol)]
#     plt.scatter(tout_all/3600, x_t)
#     plt.ylabel('mrna_'+str(gene_symbol)+'_mpc')
#     plt.xlabel('time(h)')
#     plt.ylim(0, max(x_t)*1.25)
#     plt.show

mrna_CC = list(genes_all[[12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]])

mrna_CC_data = pd.read_csv('input_files/OmicsData.txt', sep='\t', index_col=0, header=0)['Exp RNA'][[12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29]]


x_m = xoutG_all[:,282:][:, list(genes_all).index(mrna_CC[2])]


cc_dash_mrna, axs_m = plt.subplots(6,3, sharex='col')
plt.subplots_adjust(hspace = 0.4)
k=0
for i in range(6):
    for j in range(3):
        y_val = xoutG_all[:, 282:][:, list(genes_all).index(mrna_CC[k])]
        axs_m[i,j].plot(tout_all/3600, y_val)
        axs_m[i,j].axhline(y=mrna_CC_data[k],c='red')
        axs_m[i,j].set_ylim(0,max(y_val)*1.2)
        axs_m[i,j].tick_params(axis='both', which='major', labelsize='4')
        axs_m[i,j].ticklabel_format(useOffset=False, style='plain')
        axs_m[i,j].title.set_text('mrna: '+mrna_CC[k]+' (mpc)')
        axs_m[i,j].title.set_size(5)
        if i == 5:
            axs_m[i,j].set_xlabel('time(h)', fontsize=5)
        k +=1

# axs[0,2].figure

# axs[0,1].plot(tout_all/3600,xoutS_all[:, list(species_all).index('ppERK')], 'g-')
# axs[0,2].plot(tout_all/3600,xoutS_all[:, list(species_all).index('ppERK')], 'g-')
    # %% stochastic test

# flagD = 0


# flagWr = 0
# nmxlsfile = outfile

# sys.path.insert(0, os.path.abspath(model_output_dir))
# species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
#     os.path.join(input_data_folder, 'Species.txt'), encoding='latin-1')])

# species_initializations = []
# for row in species_sheet[1:]:
#     species_initializations.append(float(row[2]))

# species_initializations = np.array(species_initializations)
# species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

# model_module = importlib.import_module(model_name)
# model = model_module.getModel()
# solver = model.getSolver()          # Create solver instance
# solver.setMaxSteps = 1e10
# model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints

# xoutS_all, xoutG_all, tout_all = RunSPARCED(
#     flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)

# if flagWr == 1:
#     columnsS = [ele for ele in model.getStateIds()]
#     columnsG = columnsS[773:914]
#     resa = [sub.replace('m_', 'ag_') for sub in columnsG]
#     resi = [sub.replace('m_', 'ig_') for sub in columnsG]
#     columnsG2 = np.concatenate((resa, resi), axis=None)
#     condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS)
#     condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
#     condsSDF = None
#     condsGDF = pd.DataFrame(data=xoutG_all, columns=columnsG2)
#     condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
#     condsGDF = None



# %%
elif flagD == 1:
    flagWr = 1
    nmxlsfile = outfile

    sys.path.insert(0, os.path.abspath(model_output_dir))
    species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
        os.path.join(input_data_folder, 'Species.txt'), encoding='latin-1')])

    species_initializations = []
    for row in species_sheet[1:]:
        species_initializations.append(float(row[2]))

    species_initializations = np.array(species_initializations)
    species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0

    model_module = importlib.import_module(model_name)
    model = model_module.getModel()
    solver = model.getSolver()          # Create solver instance
    solver.setMaxSteps = 1e10
    # np.linspace(0, 30) # set timepoints
    model.setTimepoints(np.linspace(0, ts))

    xoutS_all, xoutG_all, tout_all = RunSPARCED(
        flagD, th, species_initializations, [], Vn, Vc, model)

    if flagWr == 1:
        columnsS = [ele for ele in model.getStateIds()]
        columnsG = columnsS[773:914]
        resa = [sub.replace('m_', 'ag_') for sub in columnsG]
        resi = [sub.replace('m_', 'ig_') for sub in columnsG]
        columnsG2 = np.concatenate((resa, resi), axis=None)
        condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS)
        condsSDF.to_excel(nmxlsfile+'S_0.xlsx')
        condsSDF = None
        condsGDF = pd.DataFrame(data=xoutG_all, columns=columnsG2)
        condsGDF.to_excel(nmxlsfile+'G_0.xlsx')
        condsGDF = None

# %%

# parameter check


params_check = pd.Series(data=model_module.getModel(
).getFixedParameters(), index=model.getFixedParameterIds())

params_check = params_check.iloc[params_check.index.argsort(kind='mergesort')]
params_check = params_check[:2715]

ic_check = pd.Series(data=species_initializations, index=model.getStateIds())


# %%

params_check.to_csv('params_check.csv', sep='\t')
ic_check.to_csv('ic_check.csv', sep='\t')


# %% run a few steps manually and compare


# %% write results to dataframes

# %% write results to dataframes

columnsS = [ele for ele in model.getStateIds()]
columnsG = [str('mrna_'+gene) for gene in genes_all]
resa = [sub.replace('mrna_', 'ag_') for sub in columnsG]
resi = [sub.replace('mrna_', 'ig_') for sub in columnsG]
columnsG2 = np.concatenate((resa, resi, columnsG), axis=None)
condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS, index=tout_all)
condsSDF = condsSDF.rename_axis('tout_all')
condsSDF.to_csv(nmxlsfile+'S_0.csv', sep='\t')
condsSDF = None
condsGDF = pd.DataFrame(data=xoutG_all, columns=columnsG2, index=tout_all)
condsGDF = condsGDF.rename_axis('tout_all')
condsGDF.to_csv(nmxlsfile+'G_0.csv', sep='\t')
condsGDF = None

# %%

condsGDF = pd.read_csv(nmxlsfile+'G_0.csv', sep='\t', header=0, index_col=0)

condsSDF = pd.read_csv(nmxlsfile+'S_0.csv', sep='\t', header=0, index_col=0)

# %% run odes only


sbml_file = 'SPARCED.xml'
model_name = sbml_file[0:-4]

model_output_dir = model_name
cellNumber = 0

kGsRead = pd.read_csv(os.path.join(
    input_data_folder, 'OmicsData.txt'), header=0, index_col=0, sep="\t")
gExp_mpc = np.float64(kGsRead.values[:, 0])
mExp_mpc = np.float64(kGsRead.values[:, 1])
kGin = np.float64(kGsRead.values[:, 2])
kGac = np.float64(kGsRead.values[:, 3])
kTCleak = np.float64(kGsRead.values[:, 4])
kTCmaxs = np.float64(kGsRead.values[:, 5])
kTCd = np.float64(kGsRead.values[:, 6])
mpc2nM_Vc = (1E9/(Vc*6.023E+23))

#nmxlsfile = 'GrowthStim_det_'

th = 48
Vn = 1.75E-12
Vc = 5.25E-12
# STIMligs = [100,100.0,100.0,100.0,100.0,100.0,1721.0] # EGF, Her, HGF, PDGF, FGF, IGF, INS
STIMligs = [0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


sys.path.insert(0, os.path.abspath(model_output_dir))
species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
    os.path.join(input_data_folder, 'Species.txt'), encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)
species_initializations[np.argwhere(species_initializations <= 1e-6)] = 0.0
species_initializations[155:162] = STIMligs

model_module = importlib.import_module(model_name)
model = model_module.getModel()

# %%
mrna_stable = pd.read_csv('mrna_stable.csv', sep='\t',
                          index_col=0, squeeze=True)

for k in range(len(mrna_stable)):
    model.setFixedParameterById(mrna_stable.index[k], mrna_stable[k])


# %%

paramIds = np.array(model.getFixedParameterIds())
params = np.array(model.getFixedParameters())

param_p = re.compile('mrna_\w+')

mrnaIds = np.array(list(filter(param_p.match, paramIds)))

mrna_i = np.nonzero(np.isin(paramIds, mrnaIds))[0]

params[mrna_i] = np.dot(mExp_mpc, mpc2nM_Vc)
model.setFixedParameters(params)


# %%
solver = model.getSolver()          # Create solver instance
solver.setMaxSteps = 1e10

ts = 30  # time-step to update mRNA numbers
NSteps = int(th*3600/ts)
model.setTimepoints(np.linspace(0, th, (NSteps+1)))

#xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD,th,species_initializations,[],Vn,Vc,model)
#genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs = RunPrep(flagD,Vn,model)

xoutS_all = np.zeros(shape=(NSteps+1, len(species_initializations)))
xoutS_all[0, :] = species_initializations

model.setInitialStates(xoutS_all[0, :])

rdata = amici.runAmiciSimulation(model, solver)

qq = 0

xoutS_all[(qq+1):, :] = rdata['x'][1:, :]

columnsS = [ele for ele in model.getStateIds()]
#columnsG = columnsS[773:914]
#resa = [sub.replace('m_', 'ag_') for sub in columnsG]
#resi = [sub.replace('m_', 'ig_') for sub in columnsG]
#columnsG2 = np.concatenate((resa, resi), axis=None)
condsSDF = pd.DataFrame(data=xoutS_all, columns=columnsS, index=tout_all)
condsSDF = condsSDF.rename_axis('tout_all')
condsSDF.to_csv('new_odes_nostim.csv', sep='\t')
condsSDF = None
# condsGDF = pd.DataFrame(data=xoutG_all,columns=columnsG2, index=tout_all)
# condsGDF = condsGDF.rename_axis('tout_all')
# condsGDF.to_csv(nmxlsfile+'G_0.csv', sep='\t')
# condsGDF = None

condsSDF_read = pd.read_csv('new_odes_nostim.csv',
                            sep='\t', index_col=0, header=0)

# %%

# mrna_id = model.getStateIds()[773:]
# mrna_id = [sub.replace('m_','mrna_') for sub in mrna_id]

# mrna_stable = pd.Series(data=species_initializations[773:], index=mrna_id)

# mrna_stable.to_csv('mrna_stable.csv', sep='\t', header=False)
