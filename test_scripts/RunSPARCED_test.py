#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 01:31:34 2021

@author: arnab
"""
import libsbml
import importlib
import amici
import numpy as np
import re

from modules.SGEmodule import SGEmodule
from modules.RunPrep import RunPrep

def RunSPARCED(flagD,th,spdata,genedata,Vn,Vc,model,input_data_folder):
    ts = 30 # time-step to update mRNA numbers
    NSteps = int(th*3600/ts)
    tout_all = np.arange(0,th*3600+1,30)    
    mpc2nM_Vc = (1E9/(Vc*6.023E+23))
    
    
    genedata, mExp_mpc, GenePositionMatrix, AllGenesVec, kTCmaxs, kTCleak, kTCleak2, kGin_1, kGac_1, kTCd, TARs0, tcnas, tcnrs, tck50as, tck50rs, spIDs = RunPrep(flagD,Vn,model,input_data_folder)
    
    paramIds = np.array(model.getFixedParameterIds())
    #params = np.array(model.getFixedParameters())
    
    param_p = re.compile('mrna_\w+')
    
    mrnaIds = np.array(list(filter(param_p.match,paramIds)))
    
    mrna_i = np.nonzero(np.isin(paramIds,mrnaIds))[0]
    
    species_id = list(pd.read_csv(os.path.join(input_data_folder,'Species.txt'),header=0,index_col=0,sep="\t").index)
    
    #params_new = params
    #params_new[mrna_i] = Xm
    
    xm0 = params[mrna_i]
    
    if len(spdata)==0:
        spdata0 = pd.read_csv(os.path.join(input_data_folder,'Species.txt'),header=0,index_col=0,sep="\t")
        spdata = np.float(spdata0.values[:,1])
    xoutS_all = np.zeros(shape=(NSteps+1,len(spdata)))
    xoutS_all[0,:] = spdata # 24hr time point     
    
    if len(genedata)==0:
        genedata = genedata0
    xoutG_all = np.zeros(shape=(NSteps+1,len(genedata)))
    xoutG_all[0,:] = genedata
    
    solver = model.getSolver() # Create solver instance
    solver.setMaxSteps = 1e10
    
    for qq in range(NSteps):
        params = np.array(model.getFixedParameter())
        xm0 = params[mrna_i]
        genedata,xmN,AllGenesVec = SGEmodule(flagD,ts,xoutG_all[qq,:],xoutS_all[qq,:],xm0,Vn,Vc,kTCmaxs,kTCleak,kTCd,AllGenesVec,GenePositionMatrix,kGin_1,kGac_1,tcnas,tck50as,tcnrs,tck50rs,spIDs)
        params[mnra_i] = np.dot(xmN,mpc2nM_Vc)
        model.setFixedParameters(params)
        model.setInitialStates(xoutS_all[qq,:])
        rdata = amici.runAmiciSimulation(model, solver)  # Run simulation
        xoutS_all[qq+1,:] = rdata['x'][-1,:]
        xoutG_all[qq+1,:] = genedata
        if rdata['x'][-1,species_id.index('PARP')] < rdata['x'][-1,species_id.index('cPARP')]:
            print('Apoptosis happened')
            break
    xoutS_all = xoutS_all[~np.all(xoutS_all == 0, axis=1)]
    xoutG_all = xoutG_all[~np.all(xoutG_all == 0, axis=1)]
    tout_all = tout_all[0:len(xoutS_all)]
    
    return xoutS_all, xoutG_all, tout_all