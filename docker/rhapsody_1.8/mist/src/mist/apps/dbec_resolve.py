#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sept 27, 2016

Internal pipeline for Rhapsody Analysis

@version DBECv4.3
@author: JT

"""

from scipy.stats import poisson
import numpy as np
import pandas as pd
from collections import Counter
from rpy2 import robjects
import rpy2.robjects.pandas2ri


def si_error_correction(ValidCounters, num_bc=65536, lowdepth=4, n_uniq_reads=10):
    mi_reads = [i for j in ValidCounters.values() for i in j]
    n_uniq_mi = len(set(mi_reads))
    non_sig = [i for i in mi_reads if i > 1]
    status = 'pass'
    raw_reads = sum(mi_reads)
    raw_mi = len(mi_reads)
    raw_depth = np.mean(mi_reads)
    if len(non_sig) == 0 or np.mean(non_sig) < 4:
        # no correction if only singletons or nonsingleton depth < 4
        return {'minDepth': 0, 'status': 'low_depth', 'reads': raw_reads, 'molecules': raw_mi,
                'depth': raw_depth}
    elif 1 not in mi_reads:
        # no correction if no singletons and depth >= 4
        return {'minDepth': 0, 'status': status, 'reads': raw_reads, 'molecules': raw_mi,
                'depth': raw_depth}

    count_reads = Counter(mi_reads)

    if n_uniq_mi >= n_uniq_reads:
        doublets, highinput = getDoubletFit(ValidCounters, num_bc)
        fit_reads = Counter(sorted(mi_reads, reverse=True)[doublets:])
    else:
        # add pseudo counts at 1:100 signal
        fit_reads = {i: count_reads[i]*100 if i in count_reads else 1 for i in range(max(count_reads)+1)}

    d = pd.DataFrame([[key, value] for key, value in fit_reads.items()], columns=["reads", "freq"])

    try:
        robjects.pandas2ri.activate()
        d = robjects.conversion.py2ri(d)
        params = getparams(d)
        uniq_reads = robjects.IntVector(sorted(count_reads.keys()))
        minDepth = get_min_depth(uniq_reads, params)[0]
        adj_reads = [i for i in mi_reads if i >= minDepth]
        return {'minDepth': minDepth, 'status': status, 'reads': sum(adj_reads), 'molecules': len(adj_reads),
                'depth': np.mean(adj_reads)}
    except:
        # if fail to find two peaks, remove singletons
        return {'minDepth': 2, 'status': 'separation_failed', 'reads': sum(non_sig), 'molecules': len(non_sig),
                'depth': np.mean(non_sig)}


def getDoubletRate(n_reads, bc, usePoisson=True):
    lam = float(n_reads)/bc
    return ((1-poisson.cdf(1,lam))/(1-poisson.cdf(0,lam))) if usePoisson else lam/2


def getDoubletFit(mol, bc):
    doublets = 0
    highinput = False
    for entry in mol:
        tmp_reads = mol[entry]
        reads_len = len(tmp_reads)
        if reads_len <= 400 and bc == 65536:
            # Approx >= 1 doublet at 350-400 given 4^8 MI
            continue
        d_rate = round(getDoubletRate(reads_len, bc), 3)
        d_depth = np.percentile(tmp_reads, (1 - d_rate) * 100)
        # set highinput rate to be 0.05
        if not highinput and d_rate > 0.05: highinput = True
        doublets += sum(i > d_depth for i in tmp_reads)
    #reads = sorted(reads, reverse=True)[doublets:]
    return doublets, highinput


def getSaturated(num_bc):
    if num_bc == 6561:
        return 6557
    if num_bc == 65536:
        return 65532

robjects.r('''
    options(warn=-1)
    # R: DBEC fitting
    suppressMessages(library(scales))
    suppressMessages(library(mixdist))
    suppressMessages(library(MASS))

    getparams = function(d, n = 2){
        d = as.mixdata(d[order(d$reads),])
        dpar = mixparam(c(1,d[nrow(d)-n,'reads']),c(1,1))
        params = tryCatch(fit_nbinom(d, dpar, 10), error = function(e) { NULL })
        ifelse(!is.null(params) || (n - 1) == 0, return(params), getparams(d,n-1))
    }
    fit_nbinom = function(d, dpar, size){
      # size: number of trials for testing constraints
      fit1 = mix(d, dpar,"nbinom",constr = mixconstr(consigma='NBINOM',conmu="MFX",size=rep(size,2),fixmu=c(TRUE,FALSE)))
      params = fit1$parameters
      return(params)
    }

    get_min_depth = function(uniq_reads, params){
        for(j in 1:length(uniq_reads)){
            error = dnbinom(uniq_reads[j], size = params[1,'sigma'], mu = params[1,'mu'])
            signal = dnbinom(uniq_reads[j], size = params[2,'sigma'], mu = params[2,'mu']) 
            if(signal >= error & j > 1){
                # minimum seq depth from 2peaks
                minDepth = ifelse(uniq_reads[j-1] > max(uniq_reads)/2, 2,uniq_reads[j])
                return(minDepth)
        }
        
        }
    }

    ''')

getparams = robjects.r['getparams']
get_min_depth = robjects.r['get_min_depth']