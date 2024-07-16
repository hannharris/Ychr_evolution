#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 10:55:31 2021

@author: hannahharris
"""

#!/bin/py
from Bio import SeqIO
import os
import pandas as pd
import pickle
import random
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import itertools
from scipy.spatial import distance


#get the current working directory, which is not the same as the one the file is actually in 
#https://stackoverflow.com/questions/22282760/filenotfounderror-errno-2-no-such-file-or-directory
#cwd = os.getcwd()
#files = os.listdir(cwd)
#
#
#dir= 'Users/hannahharris/Documents/Y_evolution_project/mpra_code/'
#faFile="tester.txt"
#filename = dir + faFile
#handle = open("/Users/hannahharris/Documents/Y_evolution_project/mpra_code/promoter.test.fa", "r")
#handle.readlines()


def getLocus(fasta):
    '''
    opens fasta file and extracts records 
    returns: seqRecord objects
    '''
    faFile="/Users/hannahharris/Documents/Y_evolution_project/mpra_code/" + fasta  
    print(faFile)

    handle = open(faFile, "r")
    records=list(SeqIO.parse(handle, "fasta"))
        #close the file 
    handle.close()
    #can use to print out any attributes you are curious about
    #for record in records:
    #    print(record.id)
    return records


def infoFromDescription(seqObj): 
    '''
    takes a seqObj and returns the coordinates, chromosome, and entire ID 
    '''
    #if you do nnot have strandedness info use code below:
    info= seqObj.id
    id_no_obj= info.split('(')[0]
    coords=id_no_obj.split(':')[1]
    chrom=id_no_obj.split(':')[0]
    return coords,chrom, id_no_obj



def test_for_RE_sites(element):
    '''
    takes in a string(element)
    tests for RE sites - if site exists, then return true, if doesnt exist return false 
    '''
    #patterns RE sitees XbaI, kpnI, and SfiI
    patterns =  ['TCTAGA','GGTACC', 'GGCC[AaTtCcGg]{5}GGCC']
    
    for pattern in patterns:
        if re.search(pattern, element):
            return True
    
    return False

              
#what about sequences with alot of repeating bases?

def shuffle_seq(element):
    '''
    takes in a string and returns a shuffled string
    '''
    element = list(element)
    random.shuffle(element)
    return ''.join(element)

 
def tileElements(records, tileWindow=50, tile_total_length=100):
    '''
    takes seqrecords, iterates over them and breaks them into segments
    also checks if there are any RE sites in any particular segment 
    
    returns: dictionary of sequences (values) and dictionary of rejected elements 
    '''
    element_dict={}
    rejected={}

    for record in records: ####
        locus=record.seq
        seqLen = len(locus)
        
        #tileWindow=50 #tile_size-1 
        #tile_total_length =100
        coords,chrom,id= infoFromDescription(record)
        start,end=coords.split("-")
        end = int(end)

        #make this the length that you want it to be 
        index_1= 0
        index_2 = tile_total_length 

        start_tile_coord = int(start)
        end_tile_coord = int(start) + index_2 - 1 #does the -1 fix it?


        while index_2 <= seqLen:
            element=str(locus[index_1:index_2])
            localCoords=chrom+":"+str(start_tile_coord)+"-"+str(end_tile_coord)
    
            #make sure a particular RE sequence is NOT there then add it in 
            #XbaI TCTAGA
            if test_for_RE_sites(element):
                rejected[localCoords] = element
            else: 
                element_dict[localCoords] = element
        
            #update the coordinates
            start_tile_coord += tileWindow
            end_tile_coord += tileWindow
            index_1 += tileWindow
            index_2 += tileWindow


#in the case of a not evenly divisible situation
#        if index_2 > seqLen:
#            what_to_add = seqLen - index_1 #how many extra bases are there 
#            end_tile_coord = start_tile_coord + what_to_add - 1
#            start_tile_coord = end_tile_coord - tile_total_length
#    
#            index_2  = index_1 + what_to_add
#            index_1 = index_2 - tile_total_length
#    
#            element=str(locus[index_1:index_2])
#            localCoords=chrom+":"+str(start_tile_coord)+"-"+str(end_tile_coord)
#            if test_for_RE_sites(element):
#                rejected[localCoords] = element
#            else: 
#                element_dict[localCoords] = element

    return element_dict, rejected



def generate_shuff_elements(element_dictionary):
    '''
    takes in a dictionary of elements (key: location, value: ), creates a dictionary of shuffled sequences 
    returns shuff_dict
    '''
    shuff_dict = {}
    
    for location, seq in element_dictionary.items():
        new_location = location  + "_shuff"
        new_seq = shuffle_seq(seq)
        shuff_dict[new_location] = new_seq
    
    return shuff_dict


def compare_sequences(all_potentialElements, tileWindow = 100):
    '''
    takes in dictionary of elements, aligns each element to every other element
    returns a dictionary with the key as a merge of the two aligned sequences, and the score of their alignment
    '''
  
    pairwise_dict = {}    
    
    location = list(all_potentialElements)
    sequences = list(all_potentialElements.values())
    
    location_to_compare = location[0]
    sequence_to_compare = sequences[0]
    
    del all_potentialElements[location_to_compare]

    for index in range(1, (len(all_potentialElements) +1)): #need to figure out if its actually the right length
        for a in pairwise2.align.globalxx(sequence_to_compare, sequences[index]):
            al1,al2, score, begin, end = a
        new_id=location_to_compare + "_" + location[index]
        print(score)
        if score > tileWindow:
            print("im here")
            pairwise_dict[new_id] = score
    
    
    if len(all_potentialElements.values()) == 1:
        return pairwise_dict 
    
    pairwise_dict.update(compare_sequences(all_potentialElements, 100))
   
    return pairwise_dict


def analyze_seq_comparison(pairwise_dict, tileWindow=50):
    '''
    takes in a dictionary with sequence comparisons and assess their overlap
    *sort by values
    returns an ordered dictionary by highest value to lowest
    '''
    pass
    

def add_adaptors(element_dict, tags_per_element = 5): 
    '''
    add the tags (5 per) element and upstream elements
    only includes a barcode if its within 2 mutations from another barcode 
    returns: dictionary of locations --> sequence to order
    
    '''
    barcodes=pickle.load(open("goodBarcodes_human.p","rb"))

    upstream = "ACTGGCCGCTTCACTG"
    downstream= "AGATCGGAAGAGCGTCG"
    
    all_oligos = {}
     
    barcodes_list = list(barcodes.keys())
    used =  []
    number = 0

    for key,values in element_dict.items():
        num_tags = 1  
        #add tag 5 times

        while num_tags <= tags_per_element:
            new_key = key + "_" + str(num_tags)
            tag  = random.choice(barcodes_list)
            if len(barcodes_list) == 0:
                return "ran out of barcodes"
            new_element = upstream + values + tag + downstream
            #make sure no barcode has the site 
            if test_for_RE_sites(new_element) is False:
                barcodes_list.remove(tag) #remove the choice after it's made
                distances = []
                if used:
                    for barcode in used:
                        olig_dist = distance.hamming(list(barcode), list(tag))
                        distances.append(olig_dist)
                        
                    if min(distances) >= 0.2:
                        used.append(tag) #add to the used  list
                        number += 1 
                        print("in here", number, len(barcodes_list))
                        all_oligos[new_key] = new_element
                        num_tags += 1 
                    else: 
                        print("bad oligo dist")#could delete this part - just to check that the barcodes aren't too close to  other barcodes  
                       # print(min(distances))
                       # print("failed")
                if not used: #or min(distances) <= 0.4:
                    used.append(tag) #add to the used  list 
                    all_oligos[new_key] = new_element
                    num_tags += 1             

    return all_oligos

def final_check_RE(all_oligos):
    '''
    input is a dictionary of all_oligos
    '''
    number_oligo= 0 
    for sequences in all_oligos.values():
        if test_for_RE_sites(sequences):
            number_oligo += 1
            print(sequences)
            return("there is a site", number_oligo)
    return "no sites"


if __name__ == "__main__":
    #get elements of interest

    records=getLocus('promoters_correct_apr22.fa') 
    potentialElements, rejected = tileElements(records, 50, 200)
    #print(potentialElements)
    print("length potential elements: " + str(len(potentialElements)))
    print("length rejected elements: "+ str(len(rejected)))
    
    #get positive controls
    ##pos_ctl is CMV and EEIFa1
    pos_control_records = getLocus('pos_new_apr23.fa') 
    pos_controlElements, rejected_poscontrols = tileElements(pos_control_records, 50, 200)
    print("length pos control elements: "  + str(len(pos_controlElements)))
    
    #get negative controls
    #positive control scramble, random locus scramble(actin intron?), x and y pair scramble
    neg_control_records=getLocus('all_neg_cotrol_correct_apr23.fa') 
    neg_controlElements, rejected_negcontrols = tileElements(neg_control_records,50, 200)
    shuffElements= generate_shuff_elements(neg_controlElements)
    print("length neg control elements: " + str(len(shuffElements)))

    #add the positive control dictionary and negative control dictionary to the potential elements dictionary
    all_potentialElements = {**potentialElements, **pos_controlElements}
    all_potentialElements = {**all_potentialElements, **shuffElements}
    print("length all elements: "  +  str(len(all_potentialElements)))
#    print(all_potentialElements)
    #added: 
   # print(all_potentialElements)
    #test how similar all the sequences are to each other 
    #compare_dict  = compare_sequences(all_potentialElements, 100) #change tile window inn recursive part too
   # print(compare_dict)
   
    #sliced_d = dict(itertools.islice(all_potentialElements.items(), 1000))
    #add adaptors to everything in dictionary  
#    all_oligos = add_adaptors(all_potentialElements)
#    print(all_oligos)
#    print("length all oligos: "  + str(len(all_oligos)))
#    
#    final_check_RE(all_oligos)
#    
    #save sequences that are highly overlapping 
#    df = pd.DataFrame.from_dict(compare_dict, orient='index', columns=['values'])
#    
#    df.to_csv("compare_dict.csv")
    
#    save oligos 
#    df = pd.DataFrame.from_dict(all_oligos, orient='index', columns=['elements'])
#    
#    df.to_csv("all_oligos.csv")
##    
#    #save potential elements
#    df = pd.DataFrame.from_dict(all_potentialElements, orient='index', columns=['elements'])
#    
#    df.to_csv("elements.csv")
#    
#    #save rejected elements
#    df_rejected = pd.DataFrame.from_dict(rejected, orient='index', columns=['elements'])
#    
#    df_rejected.to_csv("rejected_elements.csv")
#






