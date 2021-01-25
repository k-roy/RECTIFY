#!/usr/bin/env python

'''
oligo-dT primed RT is subject to priming at internal transcript sites with A/G-rich sequences.  
Initial inspection of pA sites mapping to coding sequences in WT cells also revealed an abundance of pA sites
upstream of lysine/arginine-rich NLS signals, which are A/G rich (>= 80% A/ <= 40 %G).  
pA sites with 6 or more bases of downstream sequence (with no CT base allowed)
that meet these criteria need to be flagged as uncertain in mapping location.  

Real 3'UTRs occasionally contain a stretch of 4A + 2G, or 5A + 1G.  As a measure to not throw away many such signals, 
one can also keep these signals for DE analysis, but keep them flagged as uncertain of true 3' ends
'''

#import sys
import GTF_GFF_manipulation, bedgraph_computation
import imp
bedgraph_computation, GTF_GFF_manipulation = imp.reload(bedgraph_computation), imp.reload(GTF_GFF_manipulation)

#print(sys.argv)
## cluster: annotated_cluster_outfilename, , feature: feature_overlaps_outfilename
#script, genome_file, GFF_FILEPATH, plus_strand_infilename, minus_strand_infilename, annotated_pA_sites_outfilename = sys.argv

genome_file="/Volumes/SPxDrive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GFF_FILEPATH="/Volumes/SPxDrive/yeast_genome_references/saccharomyces_cerevisiae_annotations_GCR1_corrected.gff"

from time import strftime
timestamp = strftime("%Y-%m-%d_%H_hr_%M_min")

def obtain_ORF_portion_of_genome():
    infile = open(GFF_FILEPATH, 'r')
    ORF_genome = {}
    for line in infile:
        info = line.strip().split('\t')
        chromosome, SGD, sequence_type, start, end, period, strand, frame, annotation = info
        if sequence_type == 'CDS'  and chromosome != 'chrmt':
            if strand == '.':
                strand = '+'
            if strand not in ORF_genome:
                ORF_genome[strand] = {}
            if chromosome not in ORF_genome[strand]:
                ORF_genome[strand][chromosome] = set([])
            start = int(start)
            end = int(end)
            for coord in range(start,end):
                ORF_genome[strand][chromosome].add(coord)         
    return ORF_genome

def load_genome(genome):
    '''
    input: genome in fasta format
    output: dictionary with chromosome names as keys, where the value is the chromosomal sequence as a string
    '''
    genome = open(genome_file, 'r')
    genome_dict = {}
    for line in genome:
        if line[0] == '>':
            current_chromosome = line.strip()[1:]
            genome_dict[current_chromosome] = ''
        else:
            genome_dict[current_chromosome] += line.strip()
    genome.close()
    print(genome_file + " genome file loaded")
    return genome_dict

def bedgraph_to_dict(filename):
    '''
    input: bedgraph, with four tab telimited columns: chromosome, start, end, reads
    output: strain/strand specific dictionary of string chromosome dictionaries of integer coordinate dictionaries with integer read values
    '''
    infile = open(filename, 'r')
    strain_dict = {}
    for line in infile:
        info = line.strip('\n').split('\t')
        chromosome, start, end, reads = info
        if chromosome not in strain_dict:
            strain_dict[chromosome] = {}
        for idx in range(int(start), int(end)):
            strain_dict[chromosome][idx] = int(float(reads))
    infile.close()
    print(filename + " bedgraph file loaded")
    return strain_dict

try:
    a = ORF_genome
except:
    ORF_genome = obtain_ORF_portion_of_genome()

try:
    a = R64_genome
except:
    R64_genome = load_genome(genome_file) 
    
try:
    a = chrom_lengths
except:
    chrom_lengths = {}
    for chrom in R64_genome:
        chrom_length = len(R64_genome[chrom])
        chrom_lengths[chrom] = chrom_length

#sample_name_lst = 'WT_rep1', 'rrp6D_rep1', 'WT_BMA64_rep1', 'rrp6D_BMA64_rep1', 'WT_BMA64_rep2', 'WT_BMA64_rep3', 'rrp6D_BMA64_rep2', 'rrp6D_BMA64_rep3',   'WT_rep2', 'WT_rep3',  'rrp6D_rep2', 'rrp6D_rep3', 'rrp6D_air1D_air2D_BMA64_rep1', 'rrp6D_air1D_air2D_BMA64_rep2', 'rrp6D_air1D_air2D_BMA64_rep3', 'rrp6D_air1D_BMA64_rep1', 'rrp6D_air1D_BMA64_rep2', 'rrp6D_air1D_BMA64_rep3', 'rrp6D_air2D_BMA64_rep1', 'rrp6D_air2D_BMA64_rep2', 'rrp6D_air2D_BMA64_rep3', 'rrp6D_BMA64_NaCl_rep1', 'rrp6D_BMA64_NaCl_rep2', 'rrp6D_BMA64_NaCl_rep3', 'rrp6D_EtOH_rep1', 'rrp6D_EtOH_rep2', 'rrp6D_EtOH_rep3', 'rrp6D_HS_rep1', 'rrp6D_HS_rep2', 'rrp6D_HS_rep3', 'WT_BMA64_NaCl_rep1', 'WT_BMA64_NaCl_rep2', 'WT_BMA64_NaCl_rep3', 'WT_EtOH_rep1', 'WT_EtOH_rep2', 'WT_EtOH_rep3', 'WT_HS_rep1', 'WT_HS_rep2', 'WT_HS_rep3'
#DIR = '/Volumes/Samsung_T5/rrp6_air_R3874/bbmap_analysis/' 


#sample_name_lst = 'RRP6_t60, RRP6_t60_ePAP, DIS3_t60, RRP43_t60'.split(', ')
#sample_name_lst = 'WT_t0, WT_t60, WT_t60_ePAP, DIS3_RRP6_t60, DIS3_RRP6_t60_ePAP'.split(', ')
sample_name_lst = 'WT_t60, DIS3_RRP6_t60'.split(', ')

DIR = '/Volumes/Samsung_T5/R3129/analysis/' 

#DIR = '/Volumes/Samsung_T5/R0356/analysis/' 
#sample_name_lst = ['WT_BMA64_37degRT_rep1_bbmap'] # , 'WT_BMA64_37degRT_rep1_bwa', 'WT_BMA64_37degRT_rep1_STAR' ] #
#sample_name_lst="WT_BMA64_37degRT_rep1_bbmap WT_BMA64_37degRT_rep2_bbmap WT_BMA64_in_vitro_polyadenylated_rep1_bbmap WT_BMA64_in_vitro_polyadenylated_rep2_bbmap".split(' ')

coord_to_info = {}

for sample_name in sample_name_lst:
    print(sample_name, 'started at:', timestamp)
    plus_strand_infilename = DIR + sample_name + '_sorted_three_prime_ends_plus_strand.bedgraph'  # _sorted
    minus_strand_infilename = DIR + sample_name + '_sorted_three_prime_ends_minus_strand.bedgraph'  # _sorted
    annotated_pA_sites_outfilename= DIR + sample_name + '_annotated_pA_sites.txt'
    
    plus_strand = bedgraph_to_dict(plus_strand_infilename)
    minus_strand = bedgraph_to_dict(minus_strand_infilename)
    annotated_pA_sites = open(annotated_pA_sites_outfilename, 'w')
    header = 'reads, chrom, coord, strand, overlaps_ORF, sequence_type, region_type, gene_id, gene_name, upstream_50bp, downstream_50bp, downstream_seq, A19, G19, A6, G6, pA_annotation, feature_annotation'.split(', ')
    header = '\t'.join(header) + '\n'
    annotated_pA_sites.write(header)
    strand = '+'
    for chrom in plus_strand:
        for coord in sorted(list(plus_strand[chrom])):
            if coord > 51 and coord < (chrom_lengths[chrom] - 51):           
                US50 = R64_genome[chrom][coord-49:coord + 1]
                DS50 = R64_genome[chrom][coord+1:coord + 51]       
                reads = plus_strand[chrom][coord] 
                pA_coord = (chrom, coord, strand)
                if pA_coord in coord_to_info:
                    complete_info = coord_to_info[pA_coord]
                else:   
                    downstream_19 = R64_genome[chrom][coord+1:coord + 20]
                    downstream_6 = downstream_19[:6]
                    A19 = downstream_19.count('A')
                    G19 = downstream_19.count('G')
                    A6 = downstream_6.count('A')
                    G6 = downstream_6.count('G')        
                    sequence_type, region_type, gene_id, gene_name = coord_to_annotation[strand][chrom].get(coord, ('intergenic', 'intergenic', 'intergenic', 'intergenic') )
                    overlaps_ORF = False
                    if coord in ORF_genome[strand][chrom]:
                        overlaps_ORF = True
                    feature_info = [chrom, coord, strand, overlaps_ORF, sequence_type, region_type, gene_id, gene_name]
                    seq_info = [US50, DS50, downstream_19, A19, G19, A6, G6]            
                    pA_annotation = '_'.join([str(e) for e in feature_info])   
                    feature_annotation = [sequence_type, region_type, gene_id, gene_name]
                    feature_annotation =  '_'.join([str(e) for e in feature_annotation])    
                    complete_info = feature_info + seq_info + [pA_annotation] + [feature_annotation]
                    coord_to_info[pA_coord] = complete_info
                output = [reads] + complete_info
                output = '\t'.join([str(e) for e in output]) + '\n'      
                annotated_pA_sites.write(output)  
    strand = '-'
    for chrom in minus_strand:
        for coord in sorted(list(minus_strand[chrom])):
            if coord > 51 and coord < (chrom_lengths[chrom] - 51):   
                US50 = bedgraph_computation.rev_comp( R64_genome[chrom][coord:coord+50] )
                DS50 = bedgraph_computation.rev_comp( R64_genome[chrom][coord-50:coord] )
                reads = minus_strand[chrom][coord]   
                pA_coord = (chrom, coord, strand)
                if pA_coord in coord_to_info:
                    complete_info = coord_to_info[pA_coord]
                else:
                    downstream_19 = bedgraph_computation.rev_comp( R64_genome[chrom][coord-19:coord] )
                    downstream_6 = bedgraph_computation.rev_comp( R64_genome[chrom][coord-6:coord] )
                    A19 = downstream_19.count('A')
                    G19 = downstream_19.count('G')
                    A6 = downstream_6.count('A')
                    G6 = downstream_6.count('G')        
                    sequence_type, region_type, gene_id, gene_name = coord_to_annotation[strand][chrom].get(coord, ('intergenic', 'intergenic', 'intergenic', 'intergenic') )
                    overlaps_ORF = False
                    if coord in ORF_genome[strand][chrom]:
                        overlaps_ORF = True
                    feature_info = [chrom, coord, strand, overlaps_ORF, sequence_type, region_type, gene_id, gene_name]
                    seq_info = [US50, DS50, downstream_19, A19, G19, A6, G6]            
                    pA_annotation = '_'.join([str(e) for e in feature_info])   
                    feature_annotation = [sequence_type, region_type, gene_id, gene_name]
                    feature_annotation =  '_'.join([str(e) for e in feature_annotation])    
                    complete_info = feature_info + seq_info + [pA_annotation] + [feature_annotation]
                    coord_to_info[pA_coord] = complete_info
                output = [reads] + complete_info
                output = '\t'.join([str(e) for e in output]) + '\n'      
                annotated_pA_sites.write(output)
                    
    annotated_pA_sites.close()
    print(sample_name, 'finished at:', timestamp)
    
#    
#    
## sample_name_lst = # 'WT_rep1', 'rrp6D_rep1', 
#sample_name_lst = 'WT_BMA64_rep1', 'rrp6D_BMA64_rep1', 'WT_BMA64_rep2', 'WT_BMA64_rep3', 'rrp6D_BMA64_rep2', 'rrp6D_BMA64_rep3',   'WT_rep2', 'WT_rep3',  'rrp6D_rep2', 'rrp6D_rep3', 'rrp6D_air1D_air2D_BMA64_rep1', 'rrp6D_air1D_air2D_BMA64_rep2', 'rrp6D_air1D_air2D_BMA64_rep3', 'rrp6D_air1D_BMA64_rep1', 'rrp6D_air1D_BMA64_rep2', 'rrp6D_air1D_BMA64_rep3', 'rrp6D_air2D_BMA64_rep1', 'rrp6D_air2D_BMA64_rep2', 'rrp6D_air2D_BMA64_rep3', 'rrp6D_BMA64_NaCl_rep1', 'rrp6D_BMA64_NaCl_rep2', 'rrp6D_BMA64_NaCl_rep3', 'rrp6D_EtOH_rep1', 'rrp6D_EtOH_rep2', 'rrp6D_EtOH_rep3', 'rrp6D_HS_rep1', 'rrp6D_HS_rep2', 'rrp6D_HS_rep3', 'WT_BMA64_NaCl_rep1', 'WT_BMA64_NaCl_rep2', 'WT_BMA64_NaCl_rep3', 'WT_EtOH_rep1', 'WT_EtOH_rep2', 'WT_EtOH_rep3', 'WT_HS_rep1', 'WT_HS_rep2', 'WT_HS_rep3'
#DIR = '/Volumes/Samsung_T5/rrp6_air_R3874/bbmap_analysis/' 
#
#
#
#for sample_name in sample_name_lst:
#    print(sample_name, 'started at:', timestamp)
#    annotated_pA_sites_infilename= DIR + sample_name + '_annotated_pA_sites.txt'
#    annotated_pA_sites_outfilename= DIR + sample_name + '_annotated_pA_sites_with_upstream_and_downstream_seq.txt'
#    annotated_pA_sites_outfile = open(annotated_pA_sites_outfilename, 'w')  
#    header = 'reads, chrom, coord, strand, overlaps_ORF, region_type, common_gene_name, systematic_gene_name, strandedness_relative_to_gene, pA_annotation, feature_annotation, upstream_50bp, downstream_50bp'.split(', ')
#    header = '\t'.join(header) + '\n'
#    annotated_pA_sites_outfile.write(header)
#    with open(annotated_pA_sites_infilename, 'r') as infile:
#        infile.readline()
#        for line in infile:
#            reads, chrom, coord, strand, overlaps_ORF, region_type, common_gene_name, systematic_gene_name, strandedness_relative_to_gene, downstream_seq, A19, G19, A6, G6, pA_annotation, feature_annotation = line.strip().split('\t')
#            coord = int(coord)
#            if coord > 100 and coord < (chrom_lengths[chrom] - 100):
#                if strand == '+':            
#                    US50 = R64_genome[chrom][coord-49:coord + 1]
#                    DS50 = R64_genome[chrom][coord+1:coord + 51]
#                else:
#                    US50 = bedgraph_computation.rev_comp( R64_genome[chrom][coord:coord+50] )
#                    DS50 = bedgraph_computation.rev_comp( R64_genome[chrom][coord-50:coord] )
#                output = [reads, chrom, str(coord), strand, overlaps_ORF, region_type, common_gene_name, systematic_gene_name, strandedness_relative_to_gene, pA_annotation, feature_annotation, US50, DS50 ]
#                output = '\t'.join([str(e) for e in output]) + '\n'            
#                annotated_pA_sites_outfile.write(output)
#infile.close()
#annotated_pA_sites_outfile.close()          