import pandas as pd
import numpy as np
import sys
import os
import re
import random
import matplotlib.pyplot as plt

from multiprocessing import Pool
from copy import deepcopy
from glob import glob
from dateutil import parser
from datetime import timedelta


###################################################
#### GLOBALS & CONSTANTS ##########################
###################################################

# Cache to reuse frequency distributions once computed
freqDists = dict()


###################################################
#### CASE DATA LOADING AND FORMATTING FUNCTION ####
###################################################

# 2010 US Census age proportions
# https://www.census.gov/prod/cen2010/briefs/c2010br-03.pdf

freqDists = dict()


def getAgeProportions(ref="USpop_byCounty_byAgegroup_2009.csv"):
    """Loads county age proportions form file"""
    countyAges = pd.read_csv(ref)
    countyAges.fips = countyAges.fips.astype(str).str.zfill(5)
    countyAges.loc[:,'label'] = countyAges.fips + countyAges.age_group
    popSums = countyAges.groupby('fips').popsize.sum()
    popPercents = {row['label']:row['popsize']/popSums[row['fips']] for index,row in countyAges.iterrows()}
    return popPercents


ageProportions = {'p':.065,
                  's':.175,
                  'a':.365,
                  'o':.264,
                  'g':.13}

countyAgeProportions = getAgeProportions()

def approximateAgeCohorts(df,
                          splits=countyAgeProportions,
                          mode='by county'):
    """Takes a df of flat cases by county and splits by age proportions"""
    dfs = []
    if mode == 'by county':
        for key in ageProportions.keys():
            dfSub = df.copy(deep=True)
            dfSub.columns = ['%s%s' % (str(i).zfill(5),key) for i in dfSub.columns]
            dfs.append(dfSub)
        dfOut = pd.concat(dfs,
                          sort=True,
                          axis=1)
        for column in dfOut.columns:
            if column not in splits.keys():
                countyKey = column[:-1]
                print("%s not found in county age data, using whole population proportions instead" % countyKey)
                defaultProportions = {countyKey+key:value for key,value in ageProportions.items()}
                splits.update(defaultProportions)
            dfOut.loc[:,column] = dfOut[column] * splits[column]
    elif mode == 'flat':
        for key,weight in splits.items():
            dfSub = df.copy(deep=True)
            dfSub.columns = ['%s%s' % (str(i).zfill(5),key) for i in dfSub.columns]
            dfSub = dfSub * weight
            dfs.append(dfSub)
        dfOut = pd.concat(dfs,
                          sort=True,
                          axis=1)
    return dfOut.sort_index()

def getSimData(scenarioRef,T=True):
    """Loads sim data"""
    df = pd.read_csv(scenarioRef,index_col=0)
    if T:
        df = df.T
    df.index = [int(i) for i in df.index]
    df.columns = [str(column).zfill(5) for column in df.columns]
    return df




########################################################
#### MODEL TRANSITION LOADING AND PARSING FUNCTIONS ####
########################################################

def dropZeros(dictionary):
    """Drops zero valued items from dictionary"""
    sansZeros = {key:value for key,value in dictionary.items() if value != 0}
    if len(sansZeros) != 0:
        dictionary = sansZeros
    return dictionary


def getRandomDistribution(row,
                    samples=50,
                    method='resample',
                    maxTries=100000,
                    minDuration=1):
    """Takes a weight, mean, and stddev, returns a probability distribution dict for temporal weighting"""
    weight = row['probability']
    mean = row['normal_mean']
    stdDev = row['normal_sd']
    random.seed('Rev. Rufus T. Cat, PhD')
    
    #Checks if said distribution has already been run
    distKey = '%s %s' % (mean,stdDev)
    
    if distKey in freqDists.keys():
        frequencies = freqDists[distKey]
    else:
        isDone = False
        tried = 0
        if method == 'resample':
            while not isDone:
                sampleSet = np.random.normal(mean,stdDev,samples)
                tried += 1
                if sampleSet.min() > minDuration-0.5:
                    isDone = True
                if tried == maxTries:
                    method = 'replace'
                    print("Exceeded max tries (%s), will now discard values < %s" % (tried,minDuration))
                    break
        if method == 'replace':
            sampleSet = []
            needed = samples
            while not isDone:
                newSamples = np.random.normal(mean,stdDev,needed)
                newSamples = [i for i in newSamples if i > minDuration]
                sampleSet += newSamples
                needed = samples-len(sampleSet)
                if needed ==0:
                    isDone = True
        sampleSet = [int(i+.5) for i in sampleSet]
        uniqueElements, countsElements = np.unique(sampleSet, return_counts=True)
        frequencies = {element:count/samples for element,count in zip(uniqueElements,countsElements)}
        freqDists[distKey] = frequencies
     
    #Weights by probability whether pulled from cache or newly generated
    frequencies = {element:weight*freq for element,freq in frequencies.items()}
    return dropZeros(frequencies)


def getDiscreteDistribution(row,discreteColumns):
    """Parses a discrete distribution table and returns a probability distribution dict for temporal weighting"""
    weight = row['probability']
    frequencies = {int(discreteColumn.split('_')[-1]):weight*row[discreteColumn] for discreteColumn in discreteColumns}
    return dropZeros(frequencies)
        

def parseTransitions(table,minDuration=1):
    """Returns the temporal weighting dict by age cohort for a given transitions table"""
    ageTransitions = dict()
    method = 'random'
    discreteColumns = [column for column in table.columns if column.startswith('discrete_')]

    if discreteColumns != []:
        method = 'discrete'

    for cohort,row in table.iterrows():
        if method == 'discrete':
            ageTransitions[cohort] = getDiscreteDistribution(row,discreteColumns)
        elif method == 'random':
            ageTransitions[cohort] = getRandomDistribution(row,minDuration=minDuration)
    return ageTransitions
            

def getTransitions(stateData,scenario,minDuration=1):
    """Loads transitions from a statedata excel file and generates a model transition representation for a given scenario"""
    print("Starting transitions prep...")
    keys = [i for i in list(pd.read_excel(stateData,None).keys()) if i.startswith('trans_')]
    flows = []
    transitions = dict()
    for key in keys:
        fromKey = key.split('_')[1]
        toKey = key.split('_')[2]
        flows.append([fromKey,toKey,key])
        table = pd.read_excel(stateData,key).fillna(0)
        if scenario not in set(table.scenario):
            lastScenario = table.scenario.tolist()[-1]
            print("Scenario %s not found for table %s, defaulting to last listed scenario instead (%s)" % (scenario,key,lastScenario))
            table = table[table.scenario==lastScenario].drop(['scenario'],axis=1).set_index('age')
        else:
            table = table[table.scenario==scenario].drop(['scenario'],axis=1).set_index('age')
        transitions[key] = parseTransitions(table,minDuration=minDuration)
    print("Done.")
    return flows,transitions


######################################################################
#### TRANSITION DISTRIBUTION OF INCIDENCE AND SUM DWELL FUNCTIONS ####
######################################################################

def dfShift(df,delay):
    """Right shifts dataframe (indexed 1 to max_day) by 'delay' days"""
    if delay == 0:
        return df
    return df.reindex(index = range(0,df.index.max()+delay+1)).shift(periods=delay,axis='index').fillna(0.0,axis=1)

def dfShiftSub(params):
    """Threaded, Right shifts dataframe (indexed 1 to max_day) by 'delay' days"""
    (df,delay,weight) = params
    if delay == 0:
        return weight*df
    return weight * df.reindex(index = range(0,df.index.max()+delay+1)).shift(periods=delay,axis='index').fillna(0.0,axis=1)


def splitAndSum(df,splits,multiThreaded):
    """Splits incidence variables by day distributions"""
    if multiThreaded:
        toProcess = [(df.copy(deep=True),delay,weight) for delay,weight in splits.items()]
        p = Pool(len(toProcess))
        results = p.map(dfShiftSub, toProcess)
        p.close()
        p.join()
    else:
        results = []
        for delay,weight in splits.items():
            splitDf = dfShift(df,delay) * weight
            results.append(splitDf)
    combinedDf = results[0]
    for splitDf in results[1:]:
        combinedDf += splitDf
    return combinedDf


def sumDwellState(dfIn,dwellTime):
    """Sums dwell times with 1 day seeing no change & accounting for fractional days"""
    dfOut = dfIn.copy(deep=True)
    while dwellTime > 1:
        if dwellTime > 2:
            increment = 1
        else:
            increment = dwellTime-1
        dfOut += dfShift(dfIn,1) * increment
        dwellTime += -1
    return dfOut

def sumDwellStateSub(params):
    """Threaded, sums dwell times with 1 day seeing no change & accounting for fractional days"""
    (dfIn,dwellTime,weight) = params

    dfOut = dfIn.copy(deep=True)
    while dwellTime > 1:
        if dwellTime > 2:
            increment = 1
        else:
            increment = dwellTime-1
        dfOut += dfShift(dfIn,1) * increment
        dwellTime += -1
    return dfOut * weight


def dwellAndSum(df,splits,multiThreaded):
    """Converts distributions of transitions out of state to cumulative dwell times"""
    if multiThreaded:
        toProcess = [(df.copy(deep=True),delay,weight) for delay,weight in splits.items()]
        p = Pool(len(toProcess))
        results = p.map(sumDwellStateSub, toProcess)
        p.close()
        p.join()
    else:
        results = []
        for delay,weight in splits.items():
            splitDf = sumDwellState(df,delay) * weight
            results.append(splitDf)
            
    combinedDf = results[0]
    for splitDf in results[1:]:
        combinedDf += splitDf
    return combinedDf



###################################################
#### OUTPUT PROCESSING FUNCTIONS ##################
###################################################


def scaleTransitions(scalers,transitionsIn):
    """Flatly adjusts transitions via dict of transition names and float scale values"""
    transitionsOut = transitionsIn
    shiftRates = lambda x: {delay:rate*scaler for delay,rate in x.items()}
    shiftAgeRates = lambda x: {agePop:shiftRates(rates) for agePop,rates in x.items()}
    for trans,scaler in scalers.items():
        transitionsOut[trans] = shiftAgeRates(transitionsOut[trans])
    return transitionsOut


def autoMerge(dfDict,toForce={'MedAttend'}):
    """Merges parallel states within an outcomes df dict by name"""
    keys = list(dfDict.keys())
    mergedDict = deepcopy(dfDict)
    for forceKey in toForce:
        for i,originalKey in enumerate(keys):
            if originalKey.endswith(forceKey):
                #print("Renaming %s to %s" % (originalKey,forceKey))
                keys[i] = forceKey
                mergedDict[forceKey] = mergedDict.pop(originalKey)
                break       
    for keyI in keys:
        for keyJ in keys:
            if keyI != keyJ:
                #print(keyI,keyJ)
                if keyJ.endswith(keyI):
                    #print("Merging %s to %s" % (keyJ,keyI))
                    mergedDict[keyI] += mergedDict[keyJ].copy(deep=True)
                    del mergedDict[keyJ]
    return mergedDict
                
        
def joinAges(dataDict):
    """Merges columns by county, dropping ages"""
    popColumns = list(dataDict.values())[0].columns.tolist()
    popColumns = [re.sub("[^0-9]", "", column) for column in popColumns]
    dictOut = dict()
    for compartmentName, table in dataDict.items():
        table.columns = popColumns
        dictOut[compartmentName] = table.sum(axis=1, level=0)
    return dictOut
    
    
def saveToCSVs(dataDict,dirOut,suffix=''):
    """Saves dict of tables to dir to csvs"""
    if not os.path.exists('data/%s' % dirOut):
        os.makedirs('data/%s' % dirOut)
    for compartmentName, table in dataDict.items():
        fileOut = 'data/%s/%s%s.csv' % (dirOut,compartmentName,suffix)
        #print("Saving table to %s" % fileOut)
        table.to_csv(fileOut)
    print('File saves done.')


###################################################
#### PRIMARY OUTCOME PROCESSING FUNCTION ##########
###################################################
    
    
def processAgeChunk(params):
    """Subprocess worker function to apply dwell time distributions"""
    (fromDf,splits,needsDwell,multiThreaded) = params
    incidenceSub = splitAndSum(fromDf,
                               splits,
                               multiThreaded)
    if needsDwell:
        dwellTotalSub = dwellAndSum(fromDf,
                                    splits,
                                    multiThreaded)
    else:
        dwellTotalSub = 'null'
    return (incidenceSub,dwellTotalSub)


def getCOVIDOutcomes(simDataRef='null',
                     scenario='null',
                     dirOut='null',
                     ignore={'trans_HypRxProt_R','trans_R_RyetS'},
                     dwellFields={'Hosp','dHosp','Vent','dVent'},
                     flowsRef='../DiseaseModel/COVID-19_CDC_Parameters.xlsx',
                     noThreading=False,
                     approximateAges=False,
                     mergeStates=True,
                     mergeAges=True,
                     saveResults=True,
                     verbose=True,
                     minDuration=1,
                     scalers='null'):
    """Takes scenario data and flow rates to generate monolithic results data structure"""
    fixBool = lambda x: str(x)[0].lower() in {'y','t'}
    
    if type(simDataRef) is str:
        if simDataRef == 'null':
            print("Error, no sim data ref passed")
            return
        simData = getSimData(simDataRef)
    else:
        simData = simDataRef.copy(deep=True)
        simDataRef = 'notebook'
        
    if scenario == 'null':
        scenario = pd.read_excel(flowsRef,'transmissibility').scenario.tolist()[-1]
        print("No scenario set, setting scenario to last in disease model (%s)" % scenario)
        
    if type(ignore) is str:
        ignore = set(ignore.split())
    
    if type(dwellFields) is str:
        dwellFields = set(dwellFields.split())
    
    approximateAges = fixBool(approximateAges)
    mergeStates = fixBool(mergeStates)
    mergeAges = fixBool(mergeAges)
    saveResults = fixBool(saveResults)
    
    if approximateAges:
        simData = approximateAgeCohorts(simData)
    
    outcomes = {'E':simData.copy(deep=True)}
    
    def addOutcome(dfIn,key):
        if key in outcomes.keys():
            outcomes[key] += dfIn
        else:
            outcomes[key] = dfIn
    
    subpops = simData.columns
    subpopAges = {subpop:subpop[-1] for subpop in subpops}
    ages = list(set(subpopAges.values()))
    allAgeSubpops = {age:[i for i in subpops if subpopAges[i] == age] for age in ages}

    flows,transitions = getTransitions(flowsRef,scenario,minDuration=minDuration)
    if scalers != 'null':
        transitions = scaleTransitions(scalers,transitions)
        
    flows = [flow for flow in flows if flow[-1] not in ignore]
    for [fromKey,toKey,tableKey] in flows:
        if verbose:
            print("From %s to %s via %s" % (fromKey,toKey,tableKey))
        fromDwellKey = fromKey + '_dwellTotal'
        incidences = []
        dwellTotals = []
        needsDwell = fromKey in dwellFields

        maxSplits = max([len(transitions[tableKey][age]) for age in ages])
        multiThreadByAges = maxSplits == 1 and not noThreading
        multiThreadChunks = (not multiThreadByAges) and (not noThreading)
        
        toProcess = [(outcomes[fromKey][allAgeSubpops[age]].copy(deep=True),
                      transitions[tableKey][age],
                      needsDwell,
                      multiThreadChunks) for age in ages]
        
        if multiThreadByAges:
            p = Pool(len(ages))
            results = p.map(processAgeChunk, toProcess)
            p.close()
            p.join()
        else:
            results = [processAgeChunk(params) for params in toProcess]

        incidences = [result[0] for result in results]
        dwellTotals = [result[1] for result in results]
            
        mergedIncidences = pd.concat(incidences,axis=1)
        mergedIncidences = mergedIncidences.reindex(subpops, axis=1)
        addOutcome(mergedIncidences,toKey)
        
        if needsDwell:
            mergedDwellTotals = pd.concat(dwellTotals,axis=1)
            mergedDwellTotals = mergedDwellTotals.reindex(subpops, axis=1)
            addOutcome(mergedDwellTotals,fromDwellKey)
    print('Outcome generation complete')
    if mergeStates:
        outcomes = autoMerge(outcomes)
    if mergeAges:
        outcomes = joinAges(outcomes)
    if saveResults:
        if dirOut == 'null':
            dirOut = simDataRef.split('/')[-1].split('.')[0]+'_'+scenario
            print("No dirOut passed, using %s" % dirOut)
        saveToCSVs(outcomes,dirOut)
    return outcomes


    

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-d', '--simDataRef', type=str, required=True,
                    help="Location of simulation date file to process")
    parser.add_argument('-s', '--scenario', type=str,
                    default='null', help='Disease model parameter set for outcome generation')
    parser.add_argument('-f', '--dirOut', type=str,
                        default = 'null',
                       help='directory under {pwd}/data/ for generated outcome file output')
    parser.add_argument('-D', '--dwellFields',
                        default = 'Hosp dHosp Vent dVent',
                       help='fields for which to cumulate dwell time, adds to compute time')
    parser.add_argument('-i', '--ignore',
                        default = 'trans_HypRxProt_R trans_R_RyetS',
                       help='loop-inducing fields to ignore from processing')
    parser.add_argument('-F', '--flowsRef', type=str,
                        default='../DiseaseModel/COVID-19_CDC_Parameters.xlsx',
                        help='location of model parameters excel file')
    parser.add_argument('-a', '--approximateAges', type=str,
                        default='y',
                        help='whether or not to split cases file by population age proportions')
    parser.add_argument('-m', '--mergeStates', type=str,
                        default='y',
                        help='whether or not to automerge similarly named compartments')
    parser.add_argument('-M', '--mergeAges', type=str,
                        default='y',
                        help='whether or not to merge counties by age cohorts')
    parser.add_argument('-m', '--minDuration', type=str,
                        default='y',
                        help='minimum transition duration for mean & std declared transitions')
   
   
    
    args = vars(parser.parse_args())
    
    getCOVIDOutcomes(**args)
