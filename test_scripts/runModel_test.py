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


species_sheet = np.array([np.array(line.strip().split("\t")) for line in open(
    os.path.join(input_data_folder, 'Species.txt'), encoding='latin-1')])

species_initializations = []
for row in species_sheet[1:]:
    species_initializations.append(float(row[2]))
species_initializations = np.array(species_initializations)

species_initializations[155:162] = STIMligs

# %% model import

sys.path.insert(0, os.path.abspath(model_output_dir))


model_module = importlib.import_module(model_name)
model = model_module.getModel()

model_param = np.array(model.getFixedParameterIds())

# %% turn off irrelevant CC reactions

ks = ['k'+str(n) for n in range(320, 432)]

ks_actual = []
for p in model_param:
    for k in ks:
        if k in p:
            ks_actual.append(p)

for p in ks_actual:
    model.setFixedParameterById(p, 0)

# %% reaction rate constant modification


model.setFixedParameterById('k424', 0.0001)  # kon, Mdi
model.setFixedParameterById('k425', 0.00001)  # koff, Mdi
# model.setFixedParameterById('k426',0.001) #kon, Mei
# model.setFixedParameterById('k427',0.0001) #koff, Mei
# model.setFixedParameterById('k428',0.001) #kon, Mai
# model.setFixedParameterById('k429',0.0001) #koff, Mai
# model.setFixedParameterById('k430',0.001) #kon, Mbi
# model.setFixedParameterById('k431',0.0001) #koff, Mbi
model.setFixedParameterById('k31_1', model.getFixedParameterById('k31_1')*2)
model.setFixedParameterById('k32_1', model.getFixedParameterById('k32_1')*2)

# model.setFixedParameterById('k12_1',model.getFixedParameterById('k12_1')/100)
# model.setFixedParameterById('k13_1',model.getFixedParameterById('k13_1')/100)
# model.setFixedParameterById('k14_1',model.getFixedParameterById('k14_1')/100)

model.setFixedParameterById('k12_1',model.getFixedParameterById('k12_1')/50)
model.setFixedParameterById('k13_1',model.getFixedParameterById('k13_1')/50)
model.setFixedParameterById('k14_1',model.getFixedParameterById('k14_1')/50)


# %%
#model_module = importlib.import_module(model_name)
#model = model_module.getModel()
flagD = 1
nmxlsfile = outfile
solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints


xoutS_all, xoutG_all, tout_all = RunSPARCED(flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)

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

# %% stochastic test
flagD = 0


flagWr = 0
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
model.setTimepoints(np.linspace(0, ts))  # np.linspace(0, 30) # set timepoints

xoutS_all, xoutG_all, tout_all = RunSPARCED(
    flagD, th, species_initializations, [], Vn, Vc, model, input_data_folder)

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

# %% test - plot trajectories
# xoutS_all[1960,:]

species_input = pd.read_csv('input_files/Species.txt', sep='\t', header=0, index_col=0)

species_all = species_input.index
genes_all = pd.read_csv('input_files/GeneReg.txt',
                        sep='\t', header=0, index_col=0).index

mpl.rcParams['figure.dpi'] = 300


def timecourse(species, x_s=xoutS_all, tout_all=tout_all):

    x_t = x_s[:, list(species_all).index(species)]
    plt.scatter(tout_all/3600, x_t)
    plt.ylabel(str(species))
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


# %%

timecourse('Cd')
timecourse('Cdk46', x_s=xoutS_all, tout_all=tout_all)
timecourse('Mdi', x_s=xoutS_all, tout_all=tout_all)
timecourse('Md', x_s=xoutS_all, tout_all=tout_all)
timecourse('ppERK')
timecourse('Mei')
timecourse('pRB')
timecourse('pRBp')
timecourse('pRBpp')
timecoruse('E2F')
# xoutS_all[:,list(species_all).index('E')]

# %%
timecourse_mrna('CCND1')
timecourse_mrna('CCND2')
timecourse_mrna('CCND3')

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
