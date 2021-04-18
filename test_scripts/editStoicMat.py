#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 13:33:51 2021

@author: arnab
"""
import pandas as pd
import numpy as np

import sys
import os

sys.path.append(os.getcwd()[0:os.getcwd().rfind('/')]+'/sparced/bin')

#%%

input_data_folder = 'input_files'


stm = pd.read_csv('input_files/StoicMat.txt', sep='\t', header=0, index_col=0)
#%%

# list all reactions that start with vC but not vCd

vC_reactions = list(filter(lambda x: x.startswith('vC') and 'vCd' not in x, stm.columns))

#vC_reactions = list(filter(lambda x: 'vC' in x and 'vCd' not in x, stm.columns))

n = list(stm.columns).index(vC_reactions[-1])

#%%

# insert zeros column at n+1

l = np.shape(stm)[0]

stm.insert(loc=n+1, column=str('vC'+str(len(vC_reactions)+1)), value=list(np.zeros(l,int)))




#%%

# save txt

stm.to_csv('input_files/StoicMat.txt', sep='\t', header=True, index=True)


#%%

# edit values

stm = pd.read_csv('input_files/StoicMat.txt', sep='\t', header=0, index_col=0)

#%%

stm.loc['Cb','vC55'] = int(1)

stm.loc['Cb_Cdk1','vC55'] = int(-1)


#%%

