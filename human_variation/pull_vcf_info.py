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

myPath = # PATH TO GITHUB

def open_genome(genome_file):
    with open(genome_file, 'r') as genome_handle:
        genome = SeqIO.to_dict(SeqIO.parse(genome_handle, 'fasta'))
        return genome
    
def extract_sequences_from_bed(bed_file, genome):
    #Function to extract sequences from a BED file
    sequences = []
    with open(bed_file, 'r') as bed_handle:
        for line in bed_handle:
            fields = line.strip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            sequence = genome[chrom][start:end]
            record = SeqRecord(sequence, id=f'{chrom}:{start}-{end}', description='')
            sequences.append(record) 
    return sequences

#from get_regions.ipynb in 1000genomes/mar24
genes = ["EIF1AX", "EIF1AY", "KDM5D" , "KDM5C","UTY", "KDM6A", "ZFY", "ZFX", "DDX3Y" ,"DDX3X", "USP9Y" , "USP9X", "RPS4Y1", "RPS4X"] 

xgenes = ["EIF1AX",  "KDM5C", "KDM6A", "ZFX", "DDX3X", "USP9X",  "RPS4X"] 

regions = ['CDS', 'intron', 'promoter'] 

for gene in genes: 
    if gene in xgenes: 
        file =  'gnomad.genomes.v4.0.sites.chrX.vcf.bgz' #/lab/solexa_page/hannah/1000genomes/gnomad/
        #https://datasetgnomad.blob.core.windows.net/dataset/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chrX.vcf.bgz
             
    else: 
        file = 'gnomad.genomes.v4.0.sites.chrY.vcf.bgz' 
        #download from: https://datasetgnomad.blob.core.windows.net/dataset/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chrY.vcf.bgz 
    
    for region in regions: 
        bedfile = myPath + "/bedfiles/" + gene + "_" + region + "_0328.bed"
        subprocess.run(['touch', "/lab/solexa_page/hannah/1000genomes/mar24/vcfs_gnomad_v4/" + gene + "." + region + ".vcf"])
        
        txt = myPath + "/vcfs_gnomad_v4/" + gene + "." + region + ".vcf"
        myfile = open(txt, 'w')
        subprocess.run(['bcftools', 'view' ,'-R', bedfile, file],  stderr=subprocess.PIPE, stdout=myfile) 
      


