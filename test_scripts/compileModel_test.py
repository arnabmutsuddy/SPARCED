#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 22:37:20 2021

@author: arnab
"""

import libsbml
import amici
import amici.plotting
import os
import numpy as np
import pandas as pd


#%%

# temporary

#input_data_folder = "input_files"

import argparse


parser = argparse.ArgumentParser(description='Provide arguments to build the SPARCED model')
parser.add_argument('--folder', metavar='folder', help='input data folder path')
parser.add_argument('--modelname', help='input sbml model name')
# parser.add_argument('--paramfile', metavar='paramfile', help='file containing any non-default values')
args = parser.parse_args()


if args.folder == None:
    print("ERROR: No input folder given. Folder must be supplied with --folder option. Use -h for more information.")
    exit(1)

if args.modelname == None:
    print("ERROR: No model name given. Model name must be supplied with --modelname option")
    exit(1)
    
if ~np.isin((str(args.modelname)+'.xml'),os.listdir()):
    print("ERROR: SBML file doesn't exit")
    exit(1)



input_data_folder = args.folder



#%% model specifications
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4] # 'BigModel_byparts_v1'
# Directory to which the generated model code is written
model_output_dir = model_name

#%% sbml

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()


#%% define observables and constant parameters
#ICf = pd.read_excel(fileSpecies,header=0,index_col=0)

ICf= pd.read_csv(os.path.join(input_data_folder,'Species.txt'),sep='\t', header = 0, index_col = 0)

ObsMat = pd.read_csv(os.path.join(input_data_folder,"Observables.txt"), sep='\t',header=0, index_col=0)
#cell_params = pd.read_csv(os.path.join(input_data_folder,"Compartments.txt"), sep='\t', header=0, index_col=0)
compartment_sheet = np.array([np.array(line.strip().split("\t")) for line in open(os.path.join(input_data_folder,'Compartments.txt'))])


#%%
Vc = float(compartment_sheet[compartment_sheet[:,0]=='Cytoplasm',1])
#Vn = cell_params.loc['Nucleus','volume']
#Vm = cell_params.loc['Mitochondrion','volume']
#Ve = cell_params.loc['Extracellular','volume']
#volumeofcell = Vc + Vn


VxPARCDL = ICf.loc[:,'compartment'].copy()
VxPARCDL = [float(compartment_sheet[compartment_sheet[:,0]==VxPARCDL[i],1]) for i in range(len(VxPARCDL))]
VxPARCDL = pd.Series(VxPARCDL, index=ICf.index)



#%%

formula_obs = []
for obs in ObsMat.columns:
    sp_obs = ObsMat.index[np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]]
    sp_obs_id = np.nonzero(np.array(ObsMat.loc[:,obs]>0))[0]
    Vr = VxPARCDL/Vc
    Vf = Vr*ObsMat.loc[:,obs].values
    
    if len(sp_obs) == 1:
        formula_i = sp_obs[0]+'*'+str(Vf[sp_obs_id][0])
    elif len(sp_obs) == 2:
        formula_i = str(sp_obs[0]+'*'+str(Vf[sp_obs_id][0])+'+'+sp_obs[1]+'*'+str(Vf[sp_obs_id][1]))
    elif len(sp_obs) > 2:
        formula_i = ''
        for j in range(len(sp_obs)-1):
            formula_i = formula_i+sp_obs[j]+'*'+str(Vf[sp_obs_id][j])+'+'
        formula_i = formula_i+str(sp_obs[-1])+'*'+str(Vf[sp_obs_id][-1])
    formula_obs.append(formula_i)

observables = {}
obs_names = list(ObsMat.columns)

for i in range(len(obs_names)):
    observables[obs_names[i]] = {}
    observables[obs_names[i]]['formula'] = formula_obs[i]





#%%


sbml_importer = amici.SbmlImporter(sbml_file)

constantParameters = [params.getId() for params in sbml_model.getListOfParameters()]

sbml_importer.sbml2amici(model_name, 
                         model_output_dir, 
                         verbose=False,
                         observables=observables,
           
                         constantParameters = constantParameters)

if np.isin(model_name,os.listdir()).all():
    print("AMICI compilation successful")
else:
    print("AMICI compilation failed")
    exit(1)
