#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 23:58:30 2021

@author: arnab
"""

import libsbml

sbml_file = 'SPARCED.xml'

sbml_reader = libsbml.SBMLReader()
sbml_doc = sbml_reader.readSBML(sbml_file)
sbml_model = sbml_doc.getModel()

#%%
sbmlns = libsbml.SBMLNamespaces()


event1 = libsbml.Event(sbmlns)

ea = libsbml.EventAssignment(sbmlns) # needs AST node math!


ea.setVariable('ppERK')
ea.setMath(libsbml.parseFormula('0'))
#libsbml.parseFormula('0.10*k4^2')

tg = libsbml.Trigger(sbmlns)
tg.setMath(libsbml.parseFormula('t - 23*3600'))
tg.setInitialValue(False)
tg.setPersistent(True)

event1.addEventAssignment(ea)
event1.setTrigger(tg)
event1.setId('e1')

sbml_model.addEvent(event1)


#%% alternative

sbmlns = libsbml.SBMLNamespaces()

event1 = sbml_model.createEvent()
event1.setUseValuesFromTriggerTime(True)

ea = event1.createEventAssignment()
ea.setVariable('ppERK')
ea.setMath(libsbml.parseFormula('0'))

tg = event1.createTrigger()
tg.setMath(libsbml.parseFormula('t - 23*3600'))
# tg.setMath(libsbml.parseFormula('t = 23*3600'))
tg.setInitialValue(False)
tg.setPersistent(True)

#%%
#required attributes?



#%%


writer = libsbml.SBMLWriter()
# writer.writeSBML(sbml_doc, str(sbml_file+'_events'))
writer.writeSBML(sbml_doc, str(sbml_file))

#%% compile with text
from antimony import *
if loadFile("SPARCED.txt") == 1:
    print("Success loading antimony file")
else:
    print("Failed to load antimony file")
    #exit(1)

if writeSBMLFile("SPARCED.xml","SPARCED") == 1:
    print("Success converting antimony to SBML")
else:
    print("Failure converting antimony to SBML")
    #exit(1)

#%%

#test events after compilation

import amici
import amici.plotting
import sys
import numpy as np
import os
import importlib
# SBML model we want to import
sbml_file = 'SPARCED.xml'
# Name of the model that will also be the name of the python module
model_name = sbml_file[0:-4]
# Directory to which the generated model code is written
model_output_dir = model_name


sys.path.insert(0, os.path.abspath(model_output_dir))


model_module = importlib.import_module(model_name)
model = model_module.getModel()

#%%

solver = model.getSolver()  # Create solver instance
solver.setMaxSteps = 1e10
model.setTimepoints(np.linspace(0, 48*3600,1000))  # np.linspace(0, 30) # set timepoints

rdata = amici.runAmiciSimulation(model, solver)
