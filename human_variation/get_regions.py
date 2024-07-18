import os
import subprocess
import Bio
from Bio.Align import MultipleSeqAlignment
import glob
from Bio import SeqRecord
from Bio.Seq import Seq #added
import pandas as pd
import genomicranges as gr
import pyranges as pr
from Bio.SeqUtils import GC
from Bio import AlignIO, SeqIO
import vcf
import random
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq 

myPath = #PATH TO GITHUB FOLDER 

def open_genome(genome_file):
    with open(genome_file, 'r') as genome_handle:
        genome = SeqIO.to_dict(SeqIO.parse(genome_handle, 'fasta'))
        print('opened genome')
        return genome

def extract_sequences_from_bed(bed_file, genome):
    
    # Function to extract sequences from a BED file
    sequences = []
    with open(bed_file, 'r') as bed_handle:
        for line in bed_handle:

            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            #print(chrom,start, end )
            sequence = genome[chrom][start:end]
            #print(sequence.seq, "TEST")
            record = SeqRecord(sequence, id=f'{chrom}:{start}-{end}', description='')
            sequences.append(record)

    return sequences

def count_num_positions(txt_name):
    #takes filepath (txt_name) and counds up the number of variants
    #returns the number of bases that have variation
    vcf_reader = vcf.Reader(filename=txt_name)

    # Iterate through the records in the VCF file
    vep_categories = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
    vep_categories = vep_categories.split('|')
    each_variant_df = pd.DataFrame(columns=['chrom', 'position', 'n_alt', 'qual', 'AC', 'AN', 'AF'] + vep_categories)

    for record in vcf_reader:
        list_of_values_for_row = [record.CHROM, record.POS, record.INFO['n_alt_alleles'], record.FILTER,  record.INFO['AC'][0], record.INFO['AN'] , record.INFO["AF"][0]] + record.INFO['vep'][0].split('|')
        each_variant_df.loc[len(each_variant_df)] = list_of_values_for_row
        
    return (each_variant_df)


all_variants = pd.DataFrame()

summary_table = pd.DataFrame(columns=['gene', 'region', 'numbervar', 'seqlen'])

genes = [ "TXLNG", "TXLNGY","EIF1AX", "EIF1AY", "KDM5D" , "KDM5C","UTY", "KDM6A", "ZFY", "ZFX", "DDX3Y" ,"DDX3X", "USP9Y" , "USP9X", "RPS4Y1", "RPS4X", "TMSB4X", "TMSB4Y", 'NLGN4X', 'NLGN4Y']
regions = ["promoter", "exon", "intron"]

for gene in genes:
    for region in regions:
        new_txt = myPath + "/vcfs_gnomad/" + gene + "." + region +  ".vcf" 
        try:
            each_variant_df = count_num_positions(new_txt)
        except:
            continue
        each_variant_df['gene'] = gene
        each_variant_df['region'] = region


        sequence_ = pd.read_table(myPath + "/bedfiles/"  + gene + "_" + region + "_0328.bed", header = None)

        sequence_.columns = ['chr', 'start', 'end']
        subtr_vals = sequence_['end']  - sequence_['start']
        sum_len = sum(subtr_vals)
        each_variant_df['sum_len'] = sum_len

        all_variants = pd.concat([all_variants, each_variant_df])

print(all_variants)

all_variants.to_csv("variants_gnomad.csv")

















