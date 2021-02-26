#### THIS FILE IS A COPY PASTE OF THE CLASSES FROM SC_VARIANTS.ipynb ON 7/16/2020
#### THIS FILE WAS LAST UPDATED ON 7/16/2020


import pandas as pd
import numpy as np
from fnmatch import fnmatch


class NetIntervention:
    def __init__(self, configs):
        self.patch_idx = None
        self.intervention_file = configs["NetInterventionFile"]

        self.src_glob = [] # array of strings
        self.dst_glob = [] # array of strings
        self.src_idxs = [] # array of array of strings
        self.dst_idxs = [] # array of array of strings
        self.T_starts = [] # array of ints
        self.T_ends = [] # array of ints
        self.fractions = [] # array of floats

        self.original_Theta = None

    def populate(self):
        with open(self.intervention_file, "r") as fobj:
            for line in fobj:
                T_start, T_end, src_glob, dst_glob, fraction = line.strip().split(" ")
                T_start, T_end = int(T_start), int(T_end)
                fraction = float(fraction)
                src_idxs = sorted(self.patch_idx[id_] for id_ in self.patch_idx if fnmatch(id_, src_glob))
                dst_idxs = sorted(self.patch_idx[id_] for id_ in self.patch_idx if fnmatch(id_, dst_glob))

                self.src_glob.append(src_glob)
                self.dst_glob.append(dst_glob)
                self.src_idxs.append(src_idxs)
                self.dst_idxs.append(dst_idxs)
                self.T_starts.append(T_start)
                self.T_ends.append(T_end)
                self.fractions.append(fraction)

    def __call__(self, configs, patch_df, params, Theta, seeds, vaxs, t, State_Array):
        if self.patch_idx is None:
            self.patch_idx = {id_: i for i, id_ in enumerate(patch_df["id"])}
            self.populate()
            self.original_Theta = Theta.copy()

        Theta[:,:,:] = self.original_Theta
        for src_idxs, dst_idxs, T_start, T_end, fraction in zip(self.src_idxs, self.dst_idxs, self.T_starts, self.T_ends, self.fractions):
            if T_start <= t < T_end:
                Theta[:, src_idxs, dst_idxs] *= fraction
                
                
class NetInterventionAdaptive:
    def __init__(self, configs):
        self.patch_idx = None
        self.intervention_file = configs["NetInterventionFile"]

        self.src_glob = [] # array of strings
        self.dst_glob = [] # array of strings
        self.src_idxs = [] # array of array of strings
        self.dst_idxs = [] # array of array of strings
        self.T_starts = [] # array of ints
        self.T_ends = [] # array of ints
        self.fractions = [] # array of floats

        self.original_Theta = None

    def populate(self):
        with open(self.intervention_file, "r") as fobj:
            for line in fobj:
                T_start, T_end, src_glob, dst_glob, fraction = line.strip().split(" ")
                T_start, T_end = int(T_start), int(T_end)
                fraction = float(fraction)
                src_idxs = sorted(self.patch_idx[id_] for id_ in self.patch_idx if fnmatch(id_, src_glob))
                dst_idxs = sorted(self.patch_idx[id_] for id_ in self.patch_idx if fnmatch(id_, dst_glob))

                self.src_glob.append(src_glob)
                self.dst_glob.append(dst_glob)
                self.src_idxs.append(src_idxs)
                self.dst_idxs.append(dst_idxs)
                self.T_starts.append(T_start)
                self.T_ends.append(T_end)
                self.fractions.append(fraction)

    def __call__(self, configs, patch_df, params, Theta, seeds, vaxs, t, State_Array):
        if self.patch_idx is None:
            self.patch_idx = {id_: i for i, id_ in enumerate(patch_df["id"])}
            self.populate()
            self.original_Theta = Theta.copy()
        
        Theta[:,:,:] = self.original_Theta
        Isymp = State_Array[-1,:,:]*float(configs['SymptomaticProbability'])
        wk_Isymp = Isymp[:-7,:].sum().sum()
        wk_Isymp = wk_Isymp/patch_df.pops.sum()
        
        for src_idxs, dst_idxs, T_start, T_end, fraction in zip(self.src_idxs, self.dst_idxs, self.T_starts, self.T_ends, self.fractions):
            if (t < T_end) & (wk_Isymp > 0.01):
                Theta[:, src_idxs, dst_idxs] *= fraction