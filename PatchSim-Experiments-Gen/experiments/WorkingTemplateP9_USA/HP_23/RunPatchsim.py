### EXPERIMENT FLAGS ITERATION

flags = {'VX2', 'EPI2', 'SC3', 'AV'}

### EXPERIMENT EXECUTION STUB

import sys
sys.path.append('../../../../PatchSim/')
sys.path.append('../../../../school_closures/')

import patchsim as sim
import pandas as pd
import numpy as np
import multiprocessing
import sc_variants as sc

from copy import deepcopy

np.random.seed(42)

### PREP SCENARIO BY FLAGS

epidemic = [flag for flag in flags if flag.startswith('EPI')][0]
schoolClosures = [flag for flag in flags if flag.startswith('SC')]

if len(schoolClosures) != 0:
    schoolClosure = schoolClosures[0]
else:
    schoolClosure = 'null'
    

    
#### SET NUMBER OF REPS AND THREADS

n = 100
threads = 50



### INITIAL LOADS OF PARAMS

print("Loading params")
configs = sim.read_config('config.patchsim')
patch_df = sim.load_patch(configs)
params = sim.load_params(configs, patch_df)
Theta = sim.load_Theta(configs, patch_df)
seeds = sim.load_seed(configs, params, patch_df)


if schoolClosure in {'SC2','SC4'}:
    scMethod = sc.NetIntervention(configs)
elif schoolClosure in {'SC1','SC3'}:
    scMethod = sc.NetInterventionAdaptive(configs)
else:
    scMethod = None


### EXPERIMENT EXECUTION

def runPatchsimSub(args):
    """Runs experiment in parallel using modified copies of params and configs"""
    (i,betaOut) = args
    print("Starting run",i)

    configsOut = deepcopy(configs)
    paramsOut = deepcopy(params)
    configsOut['ExposureRate'] = betaOut
    paramsOut['beta'] = np.where(params['beta']==1337,
                          betaOut,
                          params['beta'])
    
    df = sim.run_disease_simulation(configsOut,
                                    patch_df,
                                    params=paramsOut,
                                    Theta=Theta,
                                    seeds=seeds,
                                    write_epi=False,
                                    return_epi=True,
                                    intervene_step=scMethod)
    
    df.loc[:,'sample'] = i
    df.index.rename('id',inplace=True)
    return df


betaOut = {'EPI1':1.29e-06,
           'EPI2':1.60e-06}[epidemic]

stdDev = {'EPI1':7.08e-08,
         'EPI2':8.60e-08}[epidemic]

argsList = [(i,np.random.normal(betaOut,stdDev)) for i in range(n)]

print("Starting runs with beta %s and stddev %s" % (betaOut,stdDev))

with multiprocessing.Pool(threads) as mp_pool:
    results = mp_pool.map(runPatchsimSub, argsList)
    
results = pd.concat(results)

results.to_csv('MergedSamples.csv')