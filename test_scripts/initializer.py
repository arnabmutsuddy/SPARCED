#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import libsbml
import os
import sys
import importlib
import amici
import amici.plotting

import matplotlib.pyplot as plt
import time
import yaml
import pypesto
import pypesto.petab
import pypesto.optimize as optimize
import pypesto.visualize as visualize
import petab

#%%

# test module - load model

sbml_file = "SPARCED.xml"
model_name= sbml_file[0:-4]
model_output_dir = model_name
sys.path.insert(0, os.path.abspath(model_output_dir))
model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()
params_model = []
[params_model.append(p.getId()) for p in sbml_model.getListOfParameters()]

#%%

# get kTLids

fileRatelaws = "Ratelaws.txt"
Ratelawsf = pd.read_csv(os.path.join('input_files',fileRatelaws),header=0,index_col=0,sep='\t')
numRlaws = Ratelawsf.shape[0]

gene_params = pd.read_csv(os.path.join('input_files',"cell_definition_mcf10a_sparced.csv"), index_col="gene")
model_genes = gene_params.index
numberofgenes = len(model_genes)
S_PARCDL = pd.read_csv(os.path.join('input_files',"StoicMat.txt"), header=0, index_col=0, sep='\t')
S_TL = S_PARCDL.iloc[:,2:(2+numberofgenes)]

ObsMat = pd.read_csv(os.path.join('input_files','Observables.txt'), header=0, index_col=0, sep='\t')
NumObs = len(ObsMat.columns)


kTL_id = []
kTL_default = []
k50E_id = []
k50E_values = []

for rr in range(numRlaws):
    line1 = Ratelawsf.iloc[rr]
    if "EIF4E" in str(line1[1]):
        kTL_i = "k"+str(rr+1)+"_1"
        kTL_id.append(kTL_i)
        params_i = np.array(list(line1[2:]))
        params_i = params_i[~np.isnan(params_i)]
        kTL_value = params_i[0]
        kTL_default.append(float(kTL_value))
        a = len(params_i)       
        k50E_id_i = "k"+str(rr+1)+"_2"
        k50E_id.append(k50E_id_i)
        k50E_value_i = params_i[-1]
        k50E_values.append(k50E_value_i)        

    if Ratelawsf.index[rr] == 'vA77':
        kA77_id = "k"+str(rr+1)
    if Ratelawsf.index[rr] == 'vA87':
        kA87_id = "k"+str(rr+1)
        

kTLest = gene_params['kTLnat_s'].values
# kTLd = gene_params['kTLd_s'].values
mExp_mpc = gene_params['mrna_mpc'].copy()

cell_params = pd.read_csv(os.path.join('input_files',"Compartments.txt"), header=0, index_col=0, sep='\t')
Vc = cell_params.loc['Cytoplasm','volume']
Vn = cell_params.loc['Nucleus','volume']
Vm = cell_params.loc['Mitochondrion','volume']
Ve = cell_params.loc['Extracellular','volume']
volumeofcell = Vc + Vn


fileSpecies = 'Species.txt' # input
ICf = pd.read_csv(os.path.join('input_files',fileSpecies),header=0,index_col=0,sep='\t')
VxPARCDL = ICf.loc[:,'compartment'].copy()

for i in range(len(VxPARCDL)):
    if VxPARCDL[i] == "Cytoplasm":
        VxPARCDL[i] = Vc
    elif VxPARCDL[i] == "Nucleus":
        VxPARCDL[i] = Vn
    elif VxPARCDL[i] == "Mitochondrion":
        VxPARCDL[i] = Vm
    elif VxPARCDL[i] == "Extracellular":
        VxPARCDL[i] = Ve


VxTL = np.ones(numberofgenes)*Vc
for i in range(np.shape(S_TL)[1]):
    if len(np.nonzero(S_TL.values[:,i])[0]) != 0:
        obs_ind = int(np.nonzero(S_TL.values[:,i])[0])
        VxTL[i] = VxPARCDL[obs_ind]


xp_mpc = gene_params['prot_mpc'].copy()
pExp_nM = xp_mpc*1e9/(VxTL*6.023e23)
x0PARCDL = np.matmul(S_TL.values,pExp_nM.values)
x0PARCDL = pd.Series(x0PARCDL)
x0PARCDL.index = S_TL.index
mpc2nmcf_Vc=1E9/(Vc*6.023E+23)
mExp_nM=mExp_mpc*mpc2nmcf_Vc


x0PARCDL['Ribosome'] = ICf.IC_Xinitialized['Ribosome']
x0PARCDL['M'] = ICf.IC_Xinitialized['M']
x0PARCDL['PIP2'] = ICf.IC_Xinitialized['PIP2']
x0PARCDL['mT'] = ICf.IC_Xinitialized['mT']

k50E_default = max(k50E_values)
kA77 = 0 #BIM*Bax
kA87 = 0 #C8 activation
kC82 = 0.06/3600 ##p21 degradation
kbRi = 0
kdR0 = 0
kbR0 = 0

#%% functions
x0 = x0PARCDL

def get_observables(xout, VxPARCDL, Vc):    
    #ObsMat = pd.read_excel("observables_mat_v4.csv", header=0, index_col=0)
    Obs = []  
    Vr = VxPARCDL/Vc
    for i in range(np.shape(ObsMat)[1]):
        Obs_i = np.sum(ObsMat.values[:,i]*xout*Vr.flatten())
        Obs.append(Obs_i)
    Obs = [0 if i <= 1e-6 else i for i in Obs]
    Obs = np.array(Obs)    
    return Obs

obs0 = get_observables(x0.values, VxPARCDL.values, Vc)
def array_editor(array, inds, x0):
    for k in range(len(inds)):
        array[inds[k]] = x0[k]


def fse(x,y):    
    fs_error = ((x-y)/x)    
    for i in range(len(fs_error)):
        if np.isnan(fs_error[i]):
            fs_error[i] = 0         
    return fs_error

def fe(x,y):    
    f_error = ((x-y)/x)    
    if np.isnan(f_error):
        f_error = 0
    if x == 0 and y == 0:
        f_error = 0         
    return f_error

def timecourse(species,x_s,start_h=0,end_h=1000):    
    timeh = np.linspace(start_h,end_h,(end_h-start_h))
    species_ind = np.nonzero(S_PARCDL.index==species)[0][0]        
    x_t = x_s[start_h:end_h,species_ind]
    plt.scatter(timeh,x_t)
    plt.ylabel(str(species))
    plt.xlabel('time(h)')
    plt.title('species timecourse')
    plt.show

def timecourse_obs(obs_name,obs_all,i,start_h=0,end_h=1000,obs0_def=obs0):
    timeh = np.linspace(start_h,end_h,(end_h-start_h))
    obs_ind = np.nonzero(ObsMat.columns==str(obs_name))[0][0]
    obs_t = obs_all[start_h:end_h,obs_ind]
    plt.figure(i)
    plt.scatter(timeh,obs_t)
    plt.axhline(y=obs0_def[obs_ind],color='r')
    plt.ylabel(str(obs_name))
    plt.xlabel('time(h)')
    plt.title('obs timecourse')
    plt.savefig(os.getcwd()+'/plots/obs/'+str(obs_name)+'_'+str(time.time())+'.png', dpi = 300)
    plt.show



#%% optimizer (manual)
ts = 1000*3600
model.setTimepoints(np.linspace(0,ts,1000))
model.setFixedParameterById(kA77_id, kA77)
model.setFixedParameterById(kA87_id, kA87)

[model.setFixedParameterById(k50E_id[k],0) for k in range(len(k50E_id))]


x0 = x0PARCDL
model.setInitialStates(x0.values)
kTL_initial = kTLest

# [model.setParameterById(kTL_id[k],kTL_initial[k]) for k in range (len(kTL_id))]
[model.setFixedParameterById(kTL_id[k],kTL_initial[k]) for k in range (len(kTL_id))]


solver = model.getSolver()
solver.setMaxSteps = 1e10


rdata1 = amici.runAmiciSimulation(model,solver)
x0_1 = rdata1['x']
x1 = rdata1['x'][-1,:]
obs0_1 = rdata1['y']
obs1 = obs0_1[-1]

#adjust kTLs

error_fe = (obs0 - obs1)/obs0
kTLf_obs = np.ones(len(obs0))

for i in range(len(error_fe)):
    if error_fe[i] > 0.01 and ~np.isinf(error_fe[i]):
        kTLf_obs[i] = 1/(1-error_fe[i])
    elif error_fe[i] < -0.01 and ~np.isinf(error_fe[i]):
        kTLf_obs[i] = 1/(1-error_fe[i])
    elif error_fe[i] > -0.01 and error_fe[i] < 0.01:
        kTLf_obs[i] = 1
    # elif np.isinf(error_fe[i]):
    #     kTLf_obs[i] = 10
    elif np.isinf(error_fe[i]):
        kTLf_obs[i] = 1

kTLf = []
for i in range(len(S_TL.columns)):
    a = np.nonzero(np.array(S_TL.iloc[:,i]))[0]
    if len(a) != 0:
        sp_ind = a[0]
        obs_ind = np.nonzero(np.array(ObsMat.iloc[sp_ind,:]))[0][0]
        kTLf.append(kTLf_obs[obs_ind])
    else:
        kTLf.append(1)
kTLf = pd.Series(kTLf)
kTLf = kTLf.transform(lambda x: 1 if np.isinf(x) or np.isnan(x) else x)
kTLf = np.array(kTLf)

kTLest_new = kTLest*kTLf
kTLest_new = pd.Series(kTLest_new)
kTLest_new = kTLest_new.transform(lambda x:0 if np.isnan(x) else x)
kTLest_new = np.array(kTLest_new)
# kTLest_new[np.nonzero(np.isinf(kTLest_new))[0]]=kTLest[np.nonzero(np.isinf(kTLest_new))[0]]

# [model.setParameterById(kTL_id[k],kTLest_new[k]) for k in range(len(kTL_id))]
[model.setFixedParameterById(kTL_id[k],kTLest_new[k]) for k in range(len(kTL_id))]


rdata1_new = amici.runAmiciSimulation(model,solver)
x0_1_new = rdata1_new['x']
x1_new = rdata1_new['x'][-1]
obs0_1_new = rdata1_new['y']
obs1_new = obs0_1_new[-1]

obsfse0_1_new = fse(obs0, obs1_new)
obs_nc0_1 = ObsMat.columns[~((obsfse0_1_new > -0.01) & (obsfse0_1_new < 0.01) | (obs0==0))]

#%%

# obs_inf_nan = ObsMat.columns[np.isinf(error_fe) | np.isnan(error_fe)]

# obs_inf_i = [i for i, x in enumerate(np.isinf(error_fe) | np.isnan(error_fe)) if x]

