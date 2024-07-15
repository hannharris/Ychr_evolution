#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 14:42:23 2021

@author: hannahharris
"""
import pandas as pd
import pickle 
import feather
import os
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter


#get the working directory 
os.getcwd()
os.chdir('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results')
os.getcwd()

def read_pickle(filename):

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    return data

def read_feather(filename):

    df = feather.read_dataframe(filename)

    return df

def barcode_counts(filename):

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    return data

def get_coordinates_and_num_bcs(filtered_coords_to_barcodes):
    '''
    
    '''
    list_keys = []
    list_values = []
    num_0 = 0
    for key, value in filtered_coords_to_barcodes.items():
        print(key)
        print(value)
        value_to_add = len(value)
        print(value_to_add)

        if value_to_add == 0:
            num_0 += 1
        list_keys.append(key)
        list_values.append(value_to_add)
    
    key_value  = list(zip(list_keys,list_values))
    
    coords_num_df = pd.DataFrame(key_value, columns = ['coordinate', 'unique_bcs'], dtype = float)
    
    return coords_num_df


def overlapping_barcodes():
    '''
    '''
    pass


def make_three_df(filtered_coords_to_barcodes, barcode_counts, coords_to_bc_OG, name):
    '''
    '''
    
    summ = 0 #OG contains the raw reads
    coords_to_bc_OG_SET = {}
    for key, val in coords_to_bc_OG.items():
        summ  += len(val)
        coords_to_bc_OG_SET[key] = set(val)
        
    crs = []
    bcs = []
    number_reads = []
    
    for key, value in filtered_coords_to_barcodes.items():
        #trying to get more accurate counts 
        sample = coords_to_bc_OG[key]
        counter_s = Counter(sample)

        for bc in value: 
            crs.append(key)
            number_reads.append(counter_s[bc])
            bcs.append(bc)
            
    key_value  = list(zip(crs,bcs,number_reads))    
    coords_num_df = pd.DataFrame(key_value, columns = ['crs', 'bcs', 'number_reads'], dtype = float)
       
    name1 = name  + "_filtered1.csv"
    coords_num_df.to_csv(name1)
    
    crs = []
    bcs = []
    number_reads = []
    
    for key, value in coords_to_bc_OG_SET.items():
        sample = coords_to_bc_OG[key]
        counter_s = Counter(sample)
        
        for bc in value: 
            crs.append(key)
            number_reads.append(counter_s[bc])
            bcs.append(bc)
            
   
    key_value  = list(zip(crs,bcs,number_reads))    
    coords_num_df = pd.DataFrame(key_value, columns = ['crs', 'bcs', 'number_reads'], dtype = float)
    name2 = name + "_original1.csv"  
    coords_num_df.to_csv(name2)
    
    
    bcs = []
    number_reads = []
    for key, value, in barcode_counts.items():
        bcs.append(key)
        number_reads.append(value)
    
    key_value  = list(zip(bcs,number_reads))    
    coords_num_df = pd.DataFrame(key_value, columns = ['bcs', 'number_reads'], dtype = float)
    sum(coords_num_df['number_reads'])   
    
    name3 = name  + "_bc_reads.csv"  

    coords_num_df.to_csv(name3)
    
if __name__ == "__main__":

#    initial = read_feather("/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/bf/3a9aaef7c2dbc0cee6664dfbf20447/AGCAAGT_barcodes_per_candidate.feather")
#    no_repeats = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/bf/3a9aaef7c2dbc0cee6664dfbf20447/AGCAAGT_barcodes_per_candidate-no_repeats.feather')
#    no_repeats_noj = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/bf/3a9aaef7c2dbc0cee6664dfbf20447/AGCAAGT_barcodes_per_candidate-no_repeats-no_jackpots.feather')
    coords_to_bc_OG_AGCAAGT = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/bf/3a9aaef7c2dbc0cee6664dfbf20447/AGCAAGT_coords_to_barcodes.pickle')

#    initial = read_feather("/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/38/2ec97faaf1779604909773f51b4f99/CGGCGCT_barcodes_per_candidate.feather")
#    no_repeats = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/38/2ec97faaf1779604909773f51b4f99/CGGCGCT_barcodes_per_candidate-no_repeats.feather')
#    no_repeats_noj = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/38/2ec97faaf1779604909773f51b4f99/CGGCGCT_barcodes_per_candidate-no_repeats-no_jackpots.feather')
    coords_to_bc_OG_CGGCGCT = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/38/2ec97faaf1779604909773f51b4f99/CGGCGCT_coords_to_barcodes.pickle')
    
#    initial = read_feather("/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/60/bb9f82b2f500c747bf538a77d702aa/TCGTCAA_barcodes_per_candidate.feather")
#    no_repeats = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/60/bb9f82b2f500c747bf538a77d702aa/TCGTCAA_barcodes_per_candidate-no_repeats.feather')
#    no_repeats_noj = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/60/bb9f82b2f500c747bf538a77d702aa/TCGTCAA_barcodes_per_candidate-no_repeats-no_jackpots.feather')
    coords_to_bc_OG_TCGTCAA = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/60/bb9f82b2f500c747bf538a77d702aa/TCGTCAA_coords_to_barcodes.pickle')
    
#    initial = read_feather("/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/78/1d39d85ded8137eaa42ff24b3e16ff/AAGAAGT_barcodes_per_candidate.feather")
#    no_repeats = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/78/1d39d85ded8137eaa42ff24b3e16ff/AAGAAGT_barcodes_per_candidate-no_repeats.feather')
#    no_repeats_noj = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/78/1d39d85ded8137eaa42ff24b3e16ff/AAGAAGT_barcodes_per_candidate-no_repeats-no_jackpots.feather')
    coords_to_bc_OG_AAGAAGT = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/78/1d39d85ded8137eaa42ff24b3e16ff/AAGAAGT_coords_to_barcodes.pickle')
    
#    bf_coords_to_barcodes = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/bf/3a9aaef7c2dbc0cee6664dfbf20447/AGCAAGT_coords_to_barcodes.pickle')
#    print("here")
#    
#    initial = read_feather("/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/a9/94b51d512237526561d8c78b41c7df/TGGTCTA_barcodes_per_candidate.feather")
#    no_repeats = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/a9/94b51d512237526561d8c78b41c7df/TGGTCTA_barcodes_per_candidate-no_repeats.feather')
#    no_repeats_noj = read_feather('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/a9/94b51d512237526561d8c78b41c7df/TGGTCTA_barcodes_per_candidate-no_repeats-no_jackpots.feather')
    coords_to_bc_OG_TGGTCTA = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/a9/94b51d512237526561d8c78b41c7df/TGGTCTA_coords_to_barcodes.pickle')
 
    
    #CODE THAT DE DUPLEXES THE INFORMATION
#    summ = 0 #OG contains the raw reads
#    coords_to_bc_OG_SET = {}
#    for key, val in coords_to_bc_OG.items():
#        summ  += len(val)
#        coords_to_bc_OG_SET[key] = set(val)
#    #GETS A SET OF THE UNIQUE BCS
#    summ = 0
#    set_of_unique  = set()
#    for key, val in coords_to_bc_OG_SET.items():
#        summ += len(val)
#        for value in val:
#            set_of_unique.add(value)
#            
#    len(set_of_unique) #just the unique barcodes 
#    summ #original barcodes
#        
#    #figure out if there are overlapping BC
#    
#    sum(initial['n_barcodes'])
#    sum(no_repeats['n_barcodes']) - sum(no_repeats_noj['n_barcodes'])
#    
    #bf_coords_to_barcodes = read_pickle('/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/78/1d39d85ded8137eaa42ff24b3e16ff/AAGAAGT_coords_to_barcodes.pickle')
#    print("here")

##list_of_adaptors = ['AAGAAGTTCA', 'AGCAAGTAGT', 'CGGCGCTGGC','TCGTCAACTT', 'TGGTCTACGT']
    TGGTCTA_filtered_coords_to_barcodes = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_filtered_coords_to_barcodes.pickle")
    TGGTCTA_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_barcodes_per_candidate-no_repeats-no_jackpots.feather")
    TGGTCTA_barcode_counts = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_barcode_counts.pickle")
    
#    ####
    AAGAAGT_filtered_coords_to_barcodes = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AAGAAGT/AAGAAGT_filtered_coords_to_barcodes.pickle")
    AAGAAGT_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AAGAAGT/AAGAAGT_barcodes_per_candidate-no_repeats-no_jackpots.feather")
    AAGAAGT_barcode_counts = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AAGAAGT/AAGAAGT_barcode_counts.pickle")
######
#    
    AGCAAGT_filtered_coords_to_barcodes = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AGCAAGT/AGCAAGT_filtered_coords_to_barcodes.pickle")
    AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AGCAAGT/AGCAAGT_barcodes_per_candidate-no_repeats-no_jackpots.feather")
    AGCAAGT_barcode_counts = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/AGCAAGT/AGCAAGT_barcode_counts.pickle")
    
    CGGCGCT_filtered_coords_to_barcodes = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/CGGCGCT/CGGCGCT_filtered_coords_to_barcodes.pickle")
    CGGCGCT_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/CGGCGCT/CGGCGCT_barcodes_per_candidate-no_repeats-no_jackpots.feather")
    CGGCGCT_barcode_counts = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/CGGCGCT/CGGCGCT_barcode_counts.pickle")
    
    TCGTCAA_filtered_coords_to_barcodes = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_filtered_coords_to_barcodes.pickle")
    TCGTCAA_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_barcodes_per_candidate-no_repeats-no_jackpots.feather")
    TCGTCAA_barcode_counts = read_pickle(
            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_barcode_counts.pickle")

make_three_df(TCGTCAA_filtered_coords_to_barcodes, TCGTCAA_barcode_counts, coords_to_bc_OG_TCGTCAA, "TCGTCAACTT")
make_three_df(AGCAAGT_filtered_coords_to_barcodes, AGCAAGT_barcode_counts, coords_to_bc_OG_AGCAAGT, "AGCAAGTAGT")

make_three_df(AAGAAGT_filtered_coords_to_barcodes, AAGAAGT_barcode_counts, coords_to_bc_OG_AAGAAGT, "AAGAAGTTCA")
make_three_df(TGGTCTA_filtered_coords_to_barcodes, TGGTCTA_barcode_counts, coords_to_bc_OG_TGGTCTA, "TGGTCTACGT")
make_three_df(CGGCGCT_filtered_coords_to_barcodes, CGGCGCT_barcode_counts, coords_to_bc_OG_CGGCGCT, "CGGCGCTGGC")



#    sum(AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots['n_barcodes'])
#    
#    #filter bc based on the count to try and replicate filtereD_coords_to_barcodes (cutoff of THREE)
#    my_own_filt = {}
#    total_bc = 0
#    all_bcs = []
#    for key, val in coords_to_bc_OG_SET.items():
#        bcs_acceptable = []
#        for bc in val: 
#            if AGCAAGT_barcode_counts[bc] >= 3:
#                bcs_acceptable.append(bc)
#                all_bcs.append(AGCAAGT_barcode_counts[bc])
#        total_bc += len(bcs_acceptable)
#        my_own_filt[key] = bcs_acceptable 
#                    
#    total_bc_in_their_filtered = 0   
#    list_of_values = []
#    for key, val in AGCAAGT_filtered_coords_to_barcodes.items():
#        for value in val: 
#            list_of_values.append(AGCAAGT_barcode_counts[bc])
#        total_bc_in_their_filtered  += len(val)
#    
#    len(list_of_values)
#    max(list_of_values)
#    min(all_bcs)
#    
#    
#    #CHECK FOR OVERLAPPING BARCODES 
#    set_of_unique  = set()  
#    total = 0
#    for key, val in TCGTCAA_filtered_coords_to_barcodes.items():
#        total += len(val)
#        for value in val:
#            set_of_unique.add(value)
#    len(set_of_unique)
#    print(total)
#    #MAKE A DATAFRAME OF THE COORDS, THE BC, AND READS MAPPING TO THAT 
#    
#    
    
    
    #need to make a graph out of these **
#    sum(initial['n_barcodes'])
#    
#    summ = 0 
#    to_add = [x for x in AGCAAGT_barcode_counts.values()]
#    sum(to_add)
#    
#    summed_vals = [len(x) for x in bf_coords_to_barcodes.values()]
#    sum(summed_vals) #number of unique barcodes/things in the libary??
#    
##    summed_vals = [(x) for x in bf_coords_to_barcodes.values()]
##    sum(summed_vals) 
#    sum(AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots['n_barcodes'])
#    
#    
#    for key, value in AGCAAGT_filtered_coords_to_barcodes.items():
#        AGCAAGT_filtered_coords_to_barcodes[key] = list(value)
#    
#    
    
#    bf_coords_to_barc_NEW = {}
#    for key, values in bf_coords_to_barcodes.items():
#        bf_coords_to_barc_NEW[key] = set(values)
#    bf_coords_to_barc_NEW_counts = {}
#    for key, value in bf_coords_to_barc_NEW.items():
#        bf_coords_to_barc_NEW_counts[key] = len(value)
#        
#   sum([x for x in bf_coords_to_barc_NEW_counts.values()])
#    
#    #NOW I HAVE ALMOST THE EXACT DF AS THE AGCAAGT BARCODES PER CANDIDATE NO REPEATS NO JACKPOTS
#    AGCAAGT_coords_to_barc_NEW_counts = {}
#    for key, value in AGCAAGT_filtered_coords_to_barcodes.items():
#        AGCAAGT_coords_to_barc_NEW_counts[key] = len(value)
#   
    
#    
#    diff_dict = {}
#    diff_dict_actual_bcs = {}
#    for filt_key, filt_values in AGCAAGT_filtered_coords_to_barcodes.items():
#        for key, values in bf_coords_to_barcodes.items():  
#            if filt_key == key:
#                num_diff = len(filt_values.symmetric_difference(values))
#                diff = filt_values.symmetric_difference(values)
#                diff_dict_actual_bcs[filt_key] = diff
#                diff_dict[filt_key] = num_diff 
#      

    #make a set of all the barcodes 
    #see if the length of that is equal to the total length of barcodes 
              
#    set_bcs = set()
#    list_bcs = []
#    for key, values in bf_coords_to_barc_NEW.items():
#        for val in values:
#            set_bcs.add(val)
#            list_bcs.append(val)
#    
#    len(set_bcs)
#    len(list_bcs)
#    
#    AGCAAGT_barcode_counts['CCGAGATCAAATATG']
#    
#    #try to filter by the barcode counts/coverage and see how the numbers change 
#    try_to_match_dict = {}
#    for key, values in bf_coords_to_barc_NEW.items():
#        for value in values: 
#            if AGCAAGT_barcode_counts[value] > 50:
#                print(AGCAAGT_barcode_counts[value])
#                if key in try_to_match_dict:
#                    print('here')
#                    try_to_match_dict[key].append(value)
#                else: 
#                    try_to_match_dict[key] = [value]
#    
#    list1  =AGCAAGT_filtered_coords_to_barcodes['chr6:73520747-73520946']
#    list2 =try_to_match_dict['chr6:73520747-73520946']
#
#    len(set(list1)-set(list2))
#    len(set(list2)-set(list1))    
#    
#    union=list(set().union(list1,list2))
#    len(union)
#    
    
    #what at the files that I have and what are the outputs I need to create in order to make histograms  /visualizations of the library 
    #   file                            whats in it                                                          purpose for analysis
   # AGCAAGT_barcode_counts    counter object, each barcode with each individual read of that bc            histogram of dist. of bc reads
   # initial                    dataframe, counts of every read of each bc associated with each element     histogram od dist. of bc reads 
   # AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots  pre-filtering for bc with low counts, num unique bc per CRS                 hist. of original read counts
   # AGCAAGT_filtered_coords_to_barcodes     post-filtering for bc, num unique bc per CRS                                hist of filtered read counts 
    
#   list_keys = []
#   list_vals = []
#   for key, val in  AAGAAGT_barcode_counts.items():
#       list_keys.append(key)
#       list_vals.append(val)
#    
#   key_value  = list(zip(list_keys,list_vals))
#    
#   coords_num_df = pd.DataFrame(key_value, columns = ['coordinate', 'unique_bcs'], dtype = float)
#    
#   coords_num_df.to_csv('AAGAAGT_bc_num_df.csv')  
    
    
#   list_keys = []
#   list_vals = []
#   for key, val in  TGGTCTA_filtered_coords_to_barcodes.items():
#       list_keys.append(key)
#       list_vals.append(len(val))
#    
#   key_value  = list(zip(list_keys,list_vals))
#    
#   coords_num_df = pd.DataFrame(key_value, columns = ['coordinate', 'unique_bcs'], dtype = float)
#    
#   coords_num_df.to_csv('TGGTCTA_filtered_coords_to_barcodes_num_df.csv')
#    
#  initial.to_csv('TGGTCTA_initial.csv')  
#  TGGTCTA_barcodes_per_candidate_no_repeats_no_jackpots.to_csv('TGGTCTA_barcodes_per_candidate_no_repeats_no_jackpots_cs.csv')
#














#   AGCAAGT_barcode_counts[key]   
        
    #####      no repeats no jackpots means not by umi anymore   #####

    ####
    #what everything means 
    #AGCAAGT_barcodes_per_candidate.feather = each individual read corresponding to that barcode (~40m)
    #AGCAAGT_barcode_counts = read counts of each individual barcode (sums to the total #reads ~40m)
    #AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots = finiding the set of all barcodes (not by UMI), removing jackpotted bcs 
    
   ##AGCAAGT_filtered_coords_to_barcodes --> there must be some other filtering step to get here bc I just don't see it 

##how did they go from AGCAAGT_barcodes_per_candidate_no_repeats_no_jackpots to AGCAAGT_filtered_coords_to_barcodes
##go through it with the nums
    
   
    #try to see if there are repeated barcodes across the library ****
        #could create a set of all barcodes that would be in the feathering step ????????
        #need to get a dictionary of the barcodes per cand no rep. no filtering 

    #####  
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

#    TCGTCAA_filtered_coords_to_barcodes = read_pickle(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_filtered_coords_to_barcodes.pickle")
#    TCGTCAA_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_barcodes_per_candidate-no_repeats-no_jackpots.feather")
#    TCGTCAA_barcode_counts = read_pickle(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TCGTCAA/TCGTCAA_barcode_counts.pickle")
#    ####
#    TGGTCTA_filtered_coords_to_barcodes = read_pickle(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_filtered_coords_to_barcodes.pickle")
#    TGGTCTA_barcodes_per_candidate_no_repeats_no_jackpots = read_feather(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_barcodes_per_candidate-no_repeats-no_jackpots.feather")
#    TGGTCTA_barcode_counts = read_pickle(
#            "/Volumes/hannah/lib_mpra/mpra_assoc_sept2021/concat_results/TGGTCTA/TGGTCTA_barcode_counts.pickle")
#
#
#    CGGCGCT_coords_num_df = get_coordinates_and_num_bcs(CGGCGCT_filtered_coords_to_barcodes)
#    AAGAAGT_coords_num_df = get_coordinates_and_num_bcs(AAGAAGT_filtered_coords_to_barcodes)
#    AGCAAGT_coords_num_df = get_coordinates_and_num_bcs(AGCAAGT_filtered_coords_to_barcodes)
#    TCGTCAA_coords_num_df = get_coordinates_and_num_bcs(TCGTCAA_filtered_coords_to_barcodes)
#    TGGTCTA_coords_num_df = get_coordinates_and_num_bcs(TGGTCTA_filtered_coords_to_barcodes)
#  
#    CGGCGCT_coords_num_df.to_csv('CGGCGCT_coords_num_df.csv') 
#    AAGAAGT_coords_num_df.to_csv('AAGAAGT_coords_num_df.csv') 
#    AGCAAGT_coords_num_df.to_csv('AGCAAGT_coords_num_df.csv') 
#    TCGTCAA_coords_num_df.to_csv('TCGTCAA_coords_num_df.csv') 
#    TGGTCTA_coords_num_df.to_csv('TGGTCTA_coords_num_df.csv') 
#    
#    
#    
    

    
    

#FIGURE OUT WHICH BARCODE ARENT REPRESENTED *****
    
    


#    barcode_counter = [TGGTCTA_barcode_counts,TCGTCAA_barcode_counts,AGCAAGT_barcode_counts,
#                      AAGAAGT_barcode_counts, CGGCGCT_barcode_counts]
#    num_graphs = [1,2,3,4,5]
#    
#    graphs = []
#    for ix in range(len(barcode_counter)):
#        #ax = fig1.subplot(2,1,num_graphs[ix])
#        graph = plt_BC_hist(barcode_counter[ix]) 
#        graphs.append(graph)
#        print(graph)
#        plt.clf()
#    
#    print(graphs)
#    
#    values = []
#    for key, value in TGGTCTA_barcode_counts.items():
#        values.append(value)
#        
#    df = pd.DataFrame(values)
#    sns.histplot(df)
    
#
#AAGAAGTTCA-s_2_cat_sequence.txt.gz
#AAGAAGTTCA-s_1_cat_sequence.txt.gz
#
#AGCAAGTAGT-s_2_cat_sequence.txt.gz
#AGCAAGTAGT-s_1_cat_sequence.txt.gz
#
#CGGCGCTGGC-s_2_cat_sequence.txt.gz
#CGGCGCTGGC-s_1_cat_sequence.txt.gz
#
#TCGTCAACTT-s_2_cat_sequence.txt.gz
#TCGTCAACTT-s_1_cat_sequence.txt.gz
#
#TGGTCTACGT-s_2_cat_sequence.txt.gz
#TGGTCTACGT-s_1_cat_sequence.txt.gz
#
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AAGAAGTTCA-s_2_cat_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AAGAAGTTCA-s_1_cat_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AAGAAGTTCA-s_3_cat_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/AAGAAGTTCA
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AGCAAGTAGT-s_2_cat_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AGCAAGTAGT-s_1_cat_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/AGCAAGTAGT-s_3_cat_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/AGCAAGTAGT
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/CGGCGCTGGC-s_2_cat_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/CGGCGCTGGC-s_1_cat_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/CGGCGCTGGC-s_3_cat_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/CGGCGCTGGC
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TCGTCAACTT-s_2_cat_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TCGTCAACTT-s_1_cat_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TCGTCAACTT-s_3_cat_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/TCGTCAACTT
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_2_cat_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_1_cat_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_3_cat_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/TGGTCTACGT
#
#
#cat CGGCGCTGGC-s_2_3_sequence.txt.gz CGGCGCTGGC-s_1_3_sequence.txt.gz > CGGCGCTGGC-s_3_cat_sequence.txt.gz
#cat TCGTCAACTT-s_2_3_sequence.txt.gz TCGTCAACTT-s_1_3_sequence.txt.gz > TCGTCAACTT-s_3_cat_sequence.txt.gz
#cat TGGTCTACGT-s_2_3_sequence.txt.gz TGGTCTACGT-s_1_3_sequence.txt.gz > TGGTCTACGT-s_3_cat_sequence.txt.gz
#
##!/bin/bash
#. "/nfs/apps/test/conda_test/etc/profile.d/conda.sh"
#conda activate /lab/solexa_page/hannah/mpraflow/MPRAflow/conda-env/MPRAflow
#cd /lab/solexa_page/hannah/mpraflow/MPRAflow
#
#
#nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_1_2_sequence.txt.gz" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_1_1_sequence.txt.gz" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/TGGTCTACGT-s_1_3_sequence.txt.gz" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name assoc_basic  --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/TGGTCTACGT_one_lane

#list_of_adaptors = ['AAGAAGTTCA', 'AGCAAGTAGT', 'CGGCGCTGGC','TCGTCAACTT', 'TGGTCTACGT']
#endings= ['-s_1_1_sequence.txt.gz', '-s_2_1_sequence.txt.gz', 
#          '-s_1_2_sequence.txt.gz', '-s_2_2_sequence.txt.gz', 
#          '-s_1_3_sequence.txt.gz', '-s_2_3_sequence.txt.gz']
#commands = []
#files = []
#for adaptor in list_of_adaptors: 
#    for ending in endings:
#        pasted = "gunzip" + " " + adaptor + "" + ending
#        print(pasted)
#        ending1 = ending[:-3]
#        adaptor_ending = adaptor+ ending1
#        files.append(adaptor_ending)
#        #print(adaptor_ending)
#
#cat_files = []
#for ix in range(len(files))[::2]:
#    paste_command = "cat" + " " + files[ix] + " " + files[ix+1] + " > " + files[ix][:-12] + "cat.txt"
#    print(paste_command)
#    print("echo" ' "hi"') 
#    #gzip_command = "gzip" + " " + files[ix][:-12] + "cat.txt"
#    #print(gzip_command)
#    #to_echo ='echo "' + paste_command + '"'
#    #print(to_echo) 
#    new_file = files[ix][:-12] + "cat.txt" 
#    cat_files.append(new_file)
#    
#for ix in range(len(cat_files))[::3]:
#    file1 = cat_files[ix]
#    
#    file2 = cat_files[ix +1]
#    file3 = cat_files[ix +2]
#    print('#!/bin/bash')
#    print('. "/nfs/apps/test/conda_test/etc/profile.d/conda.sh"')
#    print('conda activate /lab/solexa_page/hannah/mpraflow/MPRAflow/conda-env/MPRAflow')
#    print('cd /lab/solexa_page/hannah/mpraflow/MPRAflow')
#    nextflow = 'nextflow run association_HH_edit.nf -w /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021 --fastq-insert "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/' + file2 + '" --fastq-insertPE "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/' + file1 + '" --fastq-bc "/lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/210909_WIGTC-HISEQ3B_HMHK3BCX3/QualityScore/' + file3 + '" --design "/lab/solexa_page/hannah/lib_mpra/lib_mpra_ordered2.fa" --name ' + cat_files[ix][:-17] + ' --outdir /lab/solexa_page/hannah/lib_mpra/mpra_assoc_sept2021/concat_results'
#    print(nextflow)
#
#    















