#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
from collections import Counter

####################
#Preparing My Data #
####################
pop = pd.read_csv("~/my.work/PhD/HomestaticExpansionProject/ModelData/NaiveTregDifferentiation.csv")
# Choosing spleen only data
Splnpop = pop.loc[pop['Organ'] == 'Spleen']

# Removing Ages that we do not have information on in the pop file
CD69Ages = [0,4,7,9,12,14,15,16,18,20]
ActivatedCD4pop = pop[pop.intage.isin(CD69Ages)].copy()

########################################
# Preparing Genevieves Activation Data #
########################################

CD69df = pd.read_csv('~/my.work/PhD/HomestaticExpansionProject/ModelData/CD69DataFromGen.csv')

# Changing IL-2-KO to KO so that the keys match properly
CD69df.loc[CD69df['Genotype'] == 'IL-2-KO', ['Genotype']] = 'KO'
#Selecting only the columns that we want
ACD69df = CD69df[[ 'Age', 'Genotype', 'CD4CD69_pct']]


#Results from the polynomial
FormResults = ([7, "WT",5.198136],
              [7, "KO", 1.642324],
              [9, "WT", 4.367956],
              [9, "KO", 3.268216])


#Grouping and finding the mean
GmCD69 = ACD69df.groupby( [ "Age", "Genotype"] ).mean().reset_index()

#Adding the results from the formula to the groupby mean results of CD69 data
for i in FormResults:
    GmCD69.loc[len(GmCD69)] = i

#########################################
# Calculating all activated CD4 T cells #
#########################################

def NumOfActivatedCD4(grp):
    '''
    This will take one pop group and finds the equivalent group from the GmCD69 group, 
    then takes the percentage from the appropriate GmCD69 percentage and multiplies it to the pop CD4CT group.
    Requirement:
         age = grp.name[0]
         Genotype = str(grp.name[1])
    '''

#The CD4Tcell groups equivalent to the CD69 data groups. This 
    pct = GmCD69.loc[(GmCD69['Age'] == grp.name[0]) & (GmCD69['Genotype']==str(grp.name[1]))]['CD4CD69_pct'].iloc[0]
    pct = pct / 100
    return grp * pct
    
# Translation of the groupby: make groups out of age and Genotype, 
#then take the CD4CT value run it through the function
ActivatedCD4pop['ActivatedCD4CT'] = ActivatedCD4pop.groupby( [ "intage", "Genotype"] )['NoTregCD4CT'].apply(NumOfActivatedCD4)

ActivatedCD4pop.to_csv('/home/jon/my.work/PhD/HomestaticExpansionProject/ModelData/ActivatedCD4pop.csv')

