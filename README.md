# STOMP-flu

This repository contains the code and resources required for STOMP-flu simulations and analyses. The included code and scenario data seeks to recreate and evaluate the effect of pharmaceutical and macroeconomic interventions upon the 2009 H1N1 influenza outbreak and hypothetical epidemics derived therefrom. 

This work was funded by US Centers for Disease Control and Prevention 75D30119C05935. 

## Model overview

STOMP-flu simulations are executed via a force of infection method within the PatchSim metapopulation patch model. This method describes the weight of exposure borne upon susceptible individuals as the sum of their daily person-hours of exposure to infectious individuals. The included scenarios contain a representation of the United States of America broken down by age and county cohorts via US Census Bureau demographic data. Exposure rates between each age by county subpopulation pair are aggregated from synthetic population activity schedule data. 

STOMP-flu's underlying disease model follows a traditional susceptible, exposed, infectious, and recovered (SEIR) structure with additional compartments to describe symptomatic treatment, hospitalization, ventilation, and death. During execution, PatchSim simulates the proportion of individuals moving through each SEIR disease state within each subpopulation. Following execution, PatchSim output data describing case counts by subpopulation, day, and iteration is analyzed via an outcome processing module to generate additional outcomes describing treatment, hospitalization, ventilation, and death.

The included scenarios explore a range of implementations of school closure, vaccination, and antiviral treatment interventions. Active school closure interventions are implemented via on the fly re-weighting of age cohort to age cohort exposure rates based upon an analysis of synthetic population activity schedules. Vaccination and antiviral interventions are implemented within the disease state model where general population and symptomatic infected individuals respectively are treated and moved into the appropriate treatment state.

The STOMP-flu project represents the work done at the Biocomplexity Institute of the University of Virginia in support of a CDC ensemble modelling effort:

[UVA STOMP-flu Project Brief](https://biocomplexity.virginia.edu/news/biocomplexity-institute-leads-cdc-research-study-influenza-mitigation-computer-modeling-and)

## Dependencies

### Libraries:
* Python 3:
* Numpy
* Pandas
* Seaborn
* Matplotlib
* Jupyterlab
* PatchSim

### Data:
* FixedUSANet.patchsim

### PatchSim configuration

In order to run the STOMP-flu experiments, the underlying PatchSim model must be either installed from github or cloned into the STOMP-flu directory.

[PatchSim Github](https://github.com/NSSAC/PatchSim)

[Tested PatchSim version on 1/27/2021](https://github.com/NSSAC/PatchSim/commit/b6d3634383e6f9f8882b52b0507d495a96a0450f)

### Data configuration:

All required data and templates required to run the STOMP-flu project are present within this repository with the exception of the underlying exposure network due to size constraints. This data file, *FixedUSANet.patchsim*, uses a synthesis of human travel models to describe exposure rates between each exposed county and age subpopulation and must be installed manually. The readme will be updated when a permanent host for this data file has been selected.

To install the population network, download and unzip *FixedUSANet.patchsim* and move it to *STOMP-flu/PatchSim-Experiments-Gen/FixedUSANet.patchsim*.


# Workflow

## Step 1: Experiment generation

**Execute PatchSim-Experiments-Gen/BuildHPExperiments.ipynb**:

This notebook parses experiment parameter files to generate ready-to-run STOMP-flu experiments: 

* Experiment scenario parameters are loaded from the passed *ExperimentSetup{experiment}.csv* file
* Flags in each scenario link to intervention code lines within *EventLines.txt*
* These flags and selected code lines are then inserted into a template experiment execution script from *templates/*
* The generated run scripts and configuration files are copied to their respective directories in *experiments/{experiment}/{scenario}/*
* Metadata of experiment parameters and file locations is written to *experiments/{experiment}/MetaData.csv*
* A monolithic bash script to run all of the experiments as well as 4 partial bash scripts for divided workflows is written to *experiments/{experiment}/RunAll.sh* and *RunAll0.sh-Runall3.sh* respectively

**Note:** at present the number of iterations and simultaneous simulation threads to run is set to 100 and 50 respectively. These may be changed as needed in by setting alternate values for the variables **n** and **threads** in the given template script.


## Step 2: Execute experiments

**Execute PatchSim-Experiments-Gen/experiments/{experiment}/RunAll.sh**

* The experiments can be run singly or by passing the 4 partial scripts to separate machines
* By default will use 50 threads and require about 5.3gb of RAM per thread
* Each iteration will require about 2 minutes per iteration with 100 iterations per scenario, 8 real life scenarios, and 24 hypothetical scenarios

**Note:** individual scenarios may be run individually by cd'ing to *Patchsim-Experiments-Gen/experiments/{experiment}/{scenario}/* and executing the python script *RunPatchSim.py.* Replicate, thread counts, and infectivity parameters may be set within this python script, while additional simulation runtime configurations are passed in *config.patchsim.* For further reference regarding PatchSim configurations, please refer to the [PatchSim documentation](https://github.com/NSSAC/PatchSim).


## Step 3: Generate outcomes

**Execute CDCFormatter/CDCFormatter.ipynb**

This notebook processes experimental exposed by subpopulation and day data to generate full outcome data across ages, geographies, and disease states:

* Disease state transitions are loaded from *FluTransitionsP3.xlsx*
* Breakdown of ages by region is loaded from *../PatchSim-Experiments-Gen/FixedUSAPop.patchsim*
* Membership of geographies by counties' first 2 FIPS digits is loaded from *SubsetRegionFIPS.csv*
* CDC output templates are loaded from *templates/{experiment-outcome}.csv*
* Experimental parameters are loaded from *../PatchSim-Experiments-Gen/experiments/{experiment}/MetaData.csv*
* The method **runAllOutcomes()** uses said metadata to load experiment scenario outputs, stochastically transition exposures between disease states, break down disease states by age and geography, and format the resulting curves in accordance to the CDC templates
* Generated data and figures are saved to the *output* and *figures* subdirectories for each scenario in *../PatchSim-Experiments-Gen/experiments/{experiment}/* 

**Note:** If **runAllOutcomes()** is run with the parameter **skipFinished=True**, then outcome processing will skip scenario subdirectories within experiments which contain an output directory. This allows for multiple instances of the notebook may be ran simultaneously without interferenceto distribute outcome processing workloads across machines.
