from GTF_GFF_manipulation import get_gene_id_from_gtf_annotation

REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI = 0.90
## check sequence for downstream AG richness and screen these out
## parameters for filtering
DOWNSTREAM_NT_TO_CHECK = 10
MAX_G_ALLOWED = 4
MAX_CT_ALLOWED = 1
## first four nucleotides must all be A to consider a read for flagging

TIFseq_filename = '/Users/kevinroy/Dropbox/temp/Pelechano_Transcript_Isoforms/S1_TIFs.txt'
#filename = '/Users/kevinroy/Dropbox/temp/Pelechano_Transcript_Isoforms/test.txt'
genome_filename = '/Users/kevinroy/Dropbox/yeast_genome/saccharomyces_cerevisiae_R64-1-1_20110208.fa'
UTR_filename = '/Users/kevinroy/Dropbox/temp/Pelechano_Transcript_Isoforms/UTR_termini_accounting_for_' + str(REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI*100) + '%_of_transcripts_screened_for_internal_AG_mispriming.txt'
gtf_filename = '/Users/kevinroy/Dropbox/yeast_genome/R64_annotations_Mito.gtf'

def load_genome(genome_filename):
    '''
    input: genome in fasta format
    output: dictionary with chromosome names as keys, where the value is the chromosomal sequence as a string
    '''
    genome = open(genome_filename, 'r')
    genome_dict = {}
    for line in genome:
        if line[0] == '>':
            current_chromosome = line.strip()[1:]
            genome_dict[current_chromosome] = ''
        else:
            genome_dict[current_chromosome] += line.strip()
    genome.close()
    print(genome_filename + " genome file loaded")
    return genome_dict

genome_dict = load_genome(genome_filename)


chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
chromosome_number_to_numeral = {}
for idx in range(16):
    chromosome_number_to_numeral[str(idx+1)] = chromosome_numerals[idx]

chromosome_lengths = {}
for chromosome in genome_dict:
    chromosome_length = len(genome_dict[chromosome])
    chromosome_lengths[chromosome] = chromosome_length

ORF_dictionary = {}
for strand in '+','-':
    ORF_dictionary[strand] = {}
flagged_isoforms = 0
total_ORF_spanning_isoforms = 0
with open(TIFseq_filename, 'r') as infile:
    for line in infile:
        line = line.strip('\n')
        if "Covering one intact ORF" in line:
            total_ORF_spanning_isoforms += 1
            idx = line.find("Y", 0, -5)
            ORF = line[idx:]
            info = line.split(' ')
            chromosome, strand, t5, t3, YPD_reads = info[:5]
            chromosome = 'chr' + chromosome_number_to_numeral[chromosome]
            chromosome_length = chromosome_lengths[chromosome]
            t5 = int(t5)
            t3 = int(t3)
            A_count, G_count, CT_count = 0,0,0
            flagged = False
            if strand == '+':
                for temp_coordinate in range(t3+1,min(chromosome_length,t3+DOWNSTREAM_NT_TO_CHECK)):
                    base = genome_dict[chromosome][temp_coordinate]
                    if base == 'A':
                        A_count += 1
                    elif base == 'G':
                        G_count += 1
                    else:
                        CT_count += 1
                if ( G_count <= MAX_G_ALLOWED and CT_count <= MAX_CT_ALLOWED ):
                    flagged = True
                    flagged_isoforms += 1
            else:
                A_count, G_count, CT_count = 0,0,0
                flagged = False
                for temp_coordinate in range(max(t3-DOWNSTREAM_NT_TO_CHECK,0),t3):
                    base = genome_dict[chromosome][temp_coordinate]
                    if base == 'T':
                        A_count += 1
                    elif base == 'C':
                        G_count += 1
                    else:
                        CT_count += 1
                if ( G_count <= MAX_G_ALLOWED and CT_count <= MAX_CT_ALLOWED ):
                    flagged = True
                    flagged_isoforms += 1
            if not flagged:
                YPD_reads = int(YPD_reads)
                if YPD_reads > 0:
                    if chromosome not in ORF_dictionary[strand]:
                        ORF_dictionary[strand][chromosome] = {}
                    if ORF not in ORF_dictionary[strand][chromosome]:
                        ORF_dictionary[strand][chromosome][ORF] = {}
                        ORF_dictionary[strand][chromosome][ORF]['5UTR'] = []
                        ORF_dictionary[strand][chromosome][ORF]['3UTR'] = []
                    ORF_dictionary[strand][chromosome][ORF]['5UTR'].append( (t5, YPD_reads))
                    ORF_dictionary[strand][chromosome][ORF]['3UTR'].append( (t3, YPD_reads))
TIFseq_filename.close()
print('ORF dictionary loaded')
print('total ORF spanning isoforms', total_ORF_spanning_isoforms)
print('total number of flagged isoforms', flagged_isoforms)

## for each ORF, order the 5'UTR and 3'UTR coordinates in order of distance from the ORF
## compute the total reads for all 5'UTR coordinates, and all 3'UTR coordinates
## starting at the ORF-proximal UTR coordinate, walk along the coordinates until REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI % of total UTR signal is accounted for
## this coordinate will be the start/end of >= REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI % transcripts from this ORF
## this value will be used in analyzing the distribution of poly(A) site data

## load the start and stop codons for each ORF
ORF_start_stop_coordinates = {}
gtf = open(gtf_filename, 'r')

ORF_start_stop_coordinates = {}
for strand in '+','-':
    ORF_start_stop_coordinates[strand] = {}
    ## load protein-coding annotations as a dictionary of strands of ORF systematic names with value as a list of lists of sequence feature types
    ## these can be of the type 3UTR, 5UTR, CDS, exon, start_codon, stop_codon followed by start and end coordinates
    ## e.g. ['5UTR', 123, 234]

for line in gtf:
    info = line.strip().split('\t')
    chromosome, region, sequence_type, start, end, period, strand, frame, annotation = info
    if 'M' not in chromosome:
        end = int(end)
        start = int(start)
        if region == 'protein_coding':
            gene_id = get_gene_id_from_gtf_annotation(annotation)
            if gene_id not in ORF_start_stop_coordinates[strand]:
                # the first entry in the ORF annotations is the chromosome
                ORF_start_stop_coordinates[strand][gene_id] = [0,0]
            if strand == '+':
                if sequence_type == 'start_codon':
                    ORF_start_stop_coordinates[strand][gene_id][0] = start
                if sequence_type == 'stop_codon':
                    ORF_start_stop_coordinates[strand][gene_id][1] = end
            elif strand == '-':
                if sequence_type == 'start_codon':
                    ORF_start_stop_coordinates[strand][gene_id][0] = end
                if sequence_type == 'stop_codon':
                    ORF_start_stop_coordinates[strand][gene_id][1] = start
        infile.close()
        print('protein_coding annotations loaded')


revised_termini = {}
for strand in '+','-':
    revised_termini[strand] = {}

    for chromosome in ORF_dictionary[strand]:
        revised_termini[strand][chromosome] = {}

        for ORF in ORF_dictionary[strand][chromosome]:
            print(ORF)
            if strand == '+':
                ORF_dictionary[strand][chromosome][ORF]['5UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['5UTR'], reverse=True)
                ORF_dictionary[strand][chromosome][ORF]['3UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['3UTR'])
            else:
                ORF_dictionary[strand][chromosome][ORF]['5UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['5UTR'])
                ORF_dictionary[strand][chromosome][ORF]['3UTR'] = sorted(ORF_dictionary[strand][chromosome][ORF]['3UTR'], reverse=True)
            fiveUTR_sum = 0
            threeUTR_sum = 0
            for coord in ORF_dictionary[strand][chromosome][ORF]['5UTR']:
                fiveUTR_sum += coord[1]
            for coord in ORF_dictionary[strand][chromosome][ORF]['3UTR']:
                threeUTR_sum += coord[1]
            fiveUTR_running_total = 0
            for coord in ORF_dictionary[strand][chromosome][ORF]['5UTR']:
                fiveUTR_running_total += coord[1]
                if fiveUTR_running_total / float(fiveUTR_sum) >= REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI:
                    fiveUTR_terminus = coord[0]
                    break
            threeUTR_running_total = 0
            for coord in ORF_dictionary[strand][chromosome][ORF]['3UTR']:
                threeUTR_running_total += coord[1]
                if threeUTR_running_total / float(threeUTR_sum) >= REQUIRED_FRACTION_FOR_TRANSCRIPT_TERMINI:
                    threeUTR_terminus = coord[0]
                    break
            revised_termini[strand][ORF] = (fiveUTR_terminus, threeUTR_terminus)

UTRs = open(UTR_filename, 'w')
# file header: strand   chrom   ORF five_UTR_start  five_UTR_end    three_UTR_start three_UTR_end
for strand in '+','-':
    for chromosome in revised_termini[strand]:
        for ORF in sorted(revised_termini[strand][chromosome]):
            ## output coordinates for five UTR start, start codon, stop codon, three UTR end
            ORF_start = ORF_start_stop_coordinates[strand][ORF][0]
            ORF_end = ORF_start_stop_coordinates[strand][ORF][1]
            UTR_info = '\t'.join([strand, chromosome,  ORF,  str(revised_termini[strand][chromosome][ORF][0]), ORF_start, ORF_end, str(revised_termini[strand][chromosome][ORF][1]) ])  + '\n'
            print(UTR_info, end=' ')
            UTRs.write(UTR_info)
UTRs.close()
