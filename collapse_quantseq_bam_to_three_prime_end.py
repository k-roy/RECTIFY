'''
Collapse quant-seq reads into a single peak to represent the three prime end
'''
import operator, pysam
print( 'psyam version: ', pysam.__version__)

## THESE FIELDS BELOW NEED TO BE ADAPTED TO USER NEEDS
genome_file='/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta'
mapping_locations="unique_and_multiple" #   unique
SAMPLE_NAMES="WT_37".split(" ") #  WT_42 WT_42_IVP
ALIGNER_SUFFIXES= "_bbmap_sorted.bam",# "_STAR_sorted_reformatted_cigar.bam", "_bwa_sorted_reformatted_cigar.bam"
ALIGNED_DIR="/Volumes/LaCie/pAmapping/aligned/"
COLLAPSED_DIR="/Volumes/LaCie/pAmapping/collapsed/"
MIN_QSCORE_FOR_START_OF_READ = 15
START_OF_READ_BP_FOR_QSCORE_FILTER = 10
MIN_BP_MATCHING_GENOME = 30
READS_TO_PROCESS = 10000 ## make this something small if troubleshooting
MIN_ALLELIC_FRACTION_FOR_POLYMORPHISM = 0.95 ## for haploid, can use high value here, diploid should be ~0.4 to account for coverage variation between two alelles
## THESE FIELDS ABOVE NEED TO BE ADAPTED TO USER NEEDS

def generate_empty_stranded_bedgraph_dict():
    '''
    takes no arguments
    returns a dictionary of strands of roman numeral-format chromosomes of empty dictionaries
    '''
    dictionary = {}
    for strand in '+-':
        dictionary[strand] = {}
        chromosome_numerals = 'I,II,III,IV,V,VI,VII,VIII,IX,X,XI,XII,XIII,XIV,XV,XVI'.split(',')
        for chromosome_numeral in chromosome_numerals:
            dictionary[strand]['chr' + chromosome_numeral] = {}
    return dictionary
    
def parse_cigar(cigar):
    '''
    input: a SAM formatted cigar string 
    output: a list of tuples, with each tuple containing the cigar operation as the first element and the number of nt as second element
    '''
    cigar_list = []
    number_of_nt_for_cigar_operation = ''
    for idx in range(len(cigar)):
        character = cigar[idx]
        if character in '0123456789':
            number_of_nt_for_cigar_operation += character
        else:
            number_of_nt = int(number_of_nt_for_cigar_operation)
            cigar_operation = character
            cigar_list.append( (cigar_operation, number_of_nt) )
            number_of_nt_for_cigar_operation = ''
    return cigar_list

def rev_comp(DNA):
    '''
    input: DNA string
    output: reverse complement string
    '''
    rev_comp_DNA = ''
    comp_bases = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', 'R':'Y', 'Y':'R', ' ':' ', 'a':'t', 't':'a', 'g':'c', 'c':'g', 'n':'n', 'r':'y', 'y':'r', '':''}
    for base in DNA[::-1]:
        rev_comp_DNA += comp_bases[base]
    return rev_comp_DNA

def write_stranded_bedgraph_dicts_to_stranded_files(bedgraph_dict, plus_strand_outfile_name, minus_strand_outfile_name, normalization_factor = 1):
    '''
    input: bedgraph dictionary, name for outfile
    output: file version of bedgraph dicitonary
    '''
    plus_strand_outfile = open(plus_strand_outfile_name, 'w')
    for chromosome in bedgraph_dict['+']:
        sorted_coordinates = sorted(bedgraph_dict['+'][chromosome].keys())
        for coordinate in sorted_coordinates:
            reads = float(bedgraph_dict['+'][chromosome][coordinate])
            #print chromosome, coordinate, reads
            output = chromosome + '\t' + str(coordinate) + '\t' + str(coordinate + 1) + '\t' + str(reads * normalization_factor) + '\n'
            #print output
            plus_strand_outfile.write(output)
    plus_strand_outfile.close()
    print(plus_strand_outfile_name + ' complete')
    minus_strand_outfile = open(minus_strand_outfile_name, 'w')
    for chromosome in bedgraph_dict['-']:
        sorted_coordinates = sorted(bedgraph_dict['-'][chromosome].keys())
        for coordinate in sorted_coordinates:
            reads = float(bedgraph_dict['-'][chromosome][coordinate])
            #print chromosome, coordinate, reads
            output = chromosome + '\t' + str(coordinate) + '\t' + str(coordinate + 1) + '\t' + str(reads * normalization_factor) + '\n'
            #print output
            minus_strand_outfile.write(output)
    minus_strand_outfile.close()
    print(minus_strand_outfile_name + ' complete')
    return None
    
def toBinary(n):
    return ''.join(str(1 & int(n) >> i) for i in range(12)[::-1])

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

try:
    a = genome_seq
except:
    genome_seq = load_genome(genome_file)

if mapping_locations == 'unique':
    MULTI_MAPPING_ALLOWED = False
elif mapping_locations == 'unique_and_multiple':
    MULTI_MAPPING_ALLOWED = True
else:
    print('you entered:', mapping_locations)
    print('you must use either "unique" or "unique_and_multiple" as the third argument')

for sample_name in SAMPLE_NAMES:
    for aligner_suffix in ALIGNER_SUFFIXES:
        infilename= ALIGNED_DIR + sample_name + aligner_suffix
        print(sample_name, aligner_suffix)
        bedgraph = generate_empty_stranded_bedgraph_dict()    
        bam_infile = pysam.AlignmentFile(infilename, "rb")
        terminal_mismatch_counts = {}
        soft_clipped_counts = {}
        three_prime_non_matching_seq_to_locations_to_counts = {}
        reads_processed = 0
        reads_collapsed = 0
        reads_removed_by_qscore = 0
        valid_alignments = 0
        reads_removed_by_min_matching_bp = 0
        flags_encountered = {}
        for read in bam_infile.fetch(until_eof=True):
            reads_processed += 1
            valid_alignment = True
            try:
                qscores = read.query_qualities
                for idx in range(START_OF_READ_BP_FOR_QSCORE_FILTER):
                    if qscores[idx] < MIN_QSCORE_FOR_START_OF_READ:
                        valid_alignment = False
                        reads_removed_by_qscore += 1
                        break
                ref_seq = read.get_reference_sequence().upper()
                flag = read.flag
                chromosome = read.reference_name
                if flag not in flags_encountered:
                    flags_encountered[flag] = 0
                flags_encountered[flag] += 1
               # NH_tag = read.get_tag("NH")
                seq = read.seq
                coordinate = read.reference_start # Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates.
                read_ID = read.query_name
                cigar = read.cigarstring 
                parsed_cigar = parse_cigar(cigar)
                mapped_matching_bp = 0
                for cigar_entry in parsed_cigar:
                    cigar_operation, number_of_nt = cigar_entry
                    if cigar_operation == '=' or cigar_operation == 'M':
                        mapped_matching_bp += number_of_nt
                if mapped_matching_bp < MIN_BP_MATCHING_GENOME:
                    valid_alignment = False
                    reads_removed_by_min_matching_bp +=1
                # if 'N' not in cigar:
                 # valid_alignment = False
            except Exception as e: 
             #   print(e)
                valid_alignment = False
#                print('here')
#                print(read)
            if  valid_alignment and 'chr' in chromosome and chromosome != 'chrMito' : # and (NH_tag == 1 or MULTI_MAPPING_ALLOWED)
                reads_collapsed += 1
                binary_flag = toBinary(flag)
                binary_flag_list = list( binary_flag )
                strand = binary_flag[-5]
                read_length = len(seq)        
                num_nt_clipped_at_end, num_nt_clipped_at_beginning  = 0, 0
                ## two possibilities for 5' discordance between read and genome
                ## 45=5S, ...AAAAG for + strand,  5S45= CTTTT.... for minus strand
                soft_clipped_or_mismatched = False
                three_prime_non_matching_seq = None
                if strand == '1':
                    strand = '+'
                    idx_in_ref = 0
                    idx_in_query = 0
                    ## if the read is mapped to the Crick strand, this read is detecting a W transcript
                    ## and the pA site is the rightmost coordinate of the mapped sequence
                    binary_flag_list[-5] = '0'  
                    in_body_of_alignment = False    
                    cigar_length = len(parsed_cigar)
                    for cigar_idx in range(cigar_length): 
                        cigar_entry = parsed_cigar[cigar_idx]
                        cigar_operation, number_of_nt = cigar_entry
                        if cigar_operation == '=' or cigar_operation == 'M':
                            in_body_of_alignment = True
                            idx_in_query += number_of_nt
                            idx_in_ref += number_of_nt
                            if cigar_idx == cigar_length-1:
                                three_prime_non_matching_seq = seq[idx_in_query-1]
                                ref_base = ref_seq[idx_in_ref-1]
                                mismatch = (three_prime_non_matching_seq, ref_base)
                                if mismatch not in terminal_mismatch_counts:
                                    terminal_mismatch_counts[mismatch] = 0
                                terminal_mismatch_counts[mismatch] += 1
                        elif cigar_operation == 'D': #  or cigar_operation == 'N':
                            idx_in_ref += number_of_nt
                        elif cigar_operation == 'I':
                            idx_in_query += number_of_nt
                        elif cigar_operation == 'X':
                            if in_body_of_alignment and number_of_nt == 1 and cigar_idx == cigar_length-1:
                                soft_clipped_or_mismatched = True
                                ## Only examine single mismatches at five prime end with an otherwise perfectly mapped read
                                three_prime_non_matching_seq = seq[idx_in_query]
                                ref_base = ref_seq[idx_in_ref]
                                mismatch = (three_prime_non_matching_seq, ref_base)
                                if mismatch not in terminal_mismatch_counts:
                                    terminal_mismatch_counts[mismatch] = 0
                                terminal_mismatch_counts[mismatch] += 1
                            idx_in_query += number_of_nt
                            idx_in_ref += number_of_nt
                        elif cigar_operation == 'S':
                            soft_clipped_or_mismatched = True
                            three_prime_non_matching_seq = seq[idx_in_query:idx_in_query+number_of_nt]
                            soft_clipped_sequence = rev_comp(three_prime_non_matching_seq)
                            if in_body_of_alignment:
                                if three_prime_non_matching_seq not in soft_clipped_counts:
                                    soft_clipped_counts[three_prime_non_matching_seq] = 0                        
                                soft_clipped_counts[three_prime_non_matching_seq] += 1
                                ## the below code block not needed as all '1S' cigar operations have been converted to '1X'
#                                if number_of_nt == 1:
#                                    three_prime_non_matching_seq = seq[idx_in_query]
#                                    ref_base = ref_seq[idx_in_ref]
#                                    mismatch = (three_prime_non_matching_seq, ref_base)
#                                    if mismatch not in terminal_mismatch_counts:
#                                        terminal_mismatch_counts[mismatch] = 0
#                                    terminal_mismatch_counts[mismatch] += 1
                            idx_in_query += number_of_nt
                        elif cigar_operation == 'N':
                            pass
                        else:
                            print('cigar operation not accounted for: ', cigar_operation, 'in', cigar, read_ID, seq)
                    adjusted_coordinate = coordinate + len(ref_seq) - 1
                    terminal_nt_in_RNA = ref_seq[idx_in_ref-1]
                    genomic_nt = genome_seq[chromosome][adjusted_coordinate]
        #            if terminal_nt_in_RNA != 'N' and terminal_nt_in_RNA != genomic_nt:
        #                print(genomic_nt,  terminal_nt_in_RNA, cigar)
        #                print(seq)
        #                print(chromosome, adjusted_coordinate)
        #                print(ref_seq, genome_seq[chromosome][adjusted_coordinate-10:adjusted_coordinate-1], genome_seq[chromosome][adjusted_coordinate-1:adjusted_coordinate+10] )        
                elif strand == '0':
                    adjusted_coordinate = coordinate
                    strand = '-'
                    ## if the read is mapped to the 0, or Watson strand, this read is detecting a Crick transcript
                    binary_flag_list[-5] = '1'
                    cigar_operation, number_of_nt = parsed_cigar[0]
                    if cigar_operation == '=':
                        three_prime_non_matching_seq = rev_comp(seq[0])
                        ref_base = rev_comp(ref_seq[0])
                        mismatch = (three_prime_non_matching_seq, ref_base)
                        if mismatch not in terminal_mismatch_counts:
                            terminal_mismatch_counts[mismatch] = 0
                        terminal_mismatch_counts[mismatch] += 1
                    elif cigar_operation == 'X':
                        if number_of_nt == 1:
                            ## Only examine single mismatches at five prime end with an otherwise perfectly mapped read
                            soft_clipped_or_mismatched = True
                            three_prime_non_matching_seq = rev_comp(seq[0])
                            ref_base = rev_comp(ref_seq[0])
                            mismatch = (three_prime_non_matching_seq, ref_base)
                            if mismatch not in terminal_mismatch_counts:
                                terminal_mismatch_counts[mismatch] = 0
                            terminal_mismatch_counts[mismatch] += 1
                    elif cigar_operation == 'S':
                        soft_clipped_or_mismatched = True
                        soft_clipped_sequence = seq[:number_of_nt]
                        three_prime_non_matching_seq = rev_comp(soft_clipped_sequence)
                        if three_prime_non_matching_seq not in soft_clipped_counts:
                            soft_clipped_counts[three_prime_non_matching_seq] = 0                        
                        soft_clipped_counts[three_prime_non_matching_seq] += 1
#                        if number_of_nt == 1:
#                            three_prime_non_matching_seq = seq[idx_in_query]
#                            ref_base = ref_seq[idx_in_ref]
#                            mismatch = (three_prime_non_matching_seq, ref_base)
#                            if mismatch not in terminal_mismatch_counts:
#                                terminal_mismatch_counts[mismatch] = 0
#                            terminal_mismatch_counts[mismatch] += 1
                    terminal_nt_in_RNA = ref_seq[0]
                    genomic_nt = genome_seq[chromosome][adjusted_coordinate ]
        #            if terminal_nt_in_RNA != genomic_nt:
        #                print(genomic_nt,  terminal_nt_in_RNA, cigar)
        #                print(seq)
        #                print(chromosome, adjusted_coordinate)
        #                print(ref_seq, genome_seq[chromosome][adjusted_coordinate-10:adjusted_coordinate-1], genome_seq[chromosome][adjusted_coordinate-1:adjusted_coordinate+10] )
                    
                    terminal_nt_in_RNA = rev_comp(terminal_nt_in_RNA)    
                    genomic_nt = rev_comp(genomic_nt)
                    
                pA_coords = (chromosome, adjusted_coordinate, strand)
#                cigar = cigar.replace('1S', '1X')      
                #   if soft_clipped_or_mismatched:  
                ## Cconvert single soft-clipped residues to single-mismatches so terminal mismatches on all aligners can be processed in the same way
                ## 1S to 1X     
                ## check for polymorphism present at allelic fraction cutoff
                if three_prime_non_matching_seq not in three_prime_non_matching_seq_to_locations_to_counts:
                    three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq] = {}
                if pA_coords not in three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq]:
                    three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq][pA_coords] = 0
                three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq][pA_coords] += 1

                if adjusted_coordinate not in bedgraph[strand][chromosome]:
                    bedgraph[strand][chromosome][adjusted_coordinate] = 0                    
                bedgraph[strand][chromosome][adjusted_coordinate] += 1    
            if valid_alignment:
                valid_alignments += 1
            if reads_processed >= READS_TO_PROCESS:
                print(reads_processed, 'reads_processed')
                print(reads_collapsed, 'reads_collapsed')
                print(valid_alignments, 'valid_alignments' )
                print(reads_removed_by_qscore, 'reads_removed_by_qscore')   
                print(reads_removed_by_min_matching_bp, 'reads_removed_by_min_matching_bp')
                print('flags_encountered:', flags_encountered)
                break
            if reads_processed % 10000 == 0:
                print(reads_removed_by_qscore, 'reads_removed_by_qscore')
                print(reads_processed, 'reads processed' )
                print(valid_alignments, 'valid_alignments' )  
                print(reads_removed_by_qscore, 'reads_removed_by_qscore')
                print('currently at ', chromosome, coordinate)    
        bam_infile.close() 
        
        
        sorted_terminal_mismatch_counts = sorted(terminal_mismatch_counts.items(), key=operator.itemgetter(1), reverse = True)
        sorted_soft_clipped_counts = sorted(soft_clipped_counts.items(), key=operator.itemgetter(1), reverse = True)
        
        filename_prefix = COLLAPSED_DIR + sample_name 
        plus_strand_outfilename = filename_prefix + '_three_prime_ends_plus_strand.bedgraph'
        minus_strand_outfilename = filename_prefix + '_three_prime_ends_minus_strand.bedgraph'
        write_stranded_bedgraph_dicts_to_stranded_files(bedgraph, plus_strand_outfilename, minus_strand_outfilename, normalization_factor = 1)
        
        soft_clipped_ends_outfile = open(filename_prefix + '_soft_clipped_ends.txt', 'w')
        soft_clipped_ends_outfile.write('soft_clipped_sequence_in_read\tsoft_clipped_sequence_in_RNA\tcounts\tfraction_of_total\ttotal_soft_clipped_counts\ttotal_uniquely_mapped_reads\ttotal_reads\n')
        total_soft_clipped_counts = 0
        for (three_prime_non_matching_seq, counts) in sorted_soft_clipped_counts:
            total_soft_clipped_counts += counts
            
        for (three_prime_non_matching_seq, counts) in sorted_soft_clipped_counts:
            fraction_of_uniquely_mapped_reads = float(counts)/reads_collapsed
            output =  [rev_comp(three_prime_non_matching_seq), three_prime_non_matching_seq,  counts, fraction_of_uniquely_mapped_reads, total_soft_clipped_counts, reads_collapsed, reads_processed]
            output = [str(e) for e in output]
            soft_clipped_ends_outfile.write('\t'.join(output) + '\n')
        soft_clipped_ends_outfile.close()
        
        single_mismatched_ends_outfile = open(filename_prefix + '_single_mismatched_ends.txt', 'w')
        single_mismatched_ends_outfile.write('last_sequenced_base_before_polyA_tail\tgenomic_base\tcounts\tfraction_total_counts\ttotal_counts\n')
        total_single_mismatched_counts = 0
        for (mismatch_pair, counts) in sorted_terminal_mismatch_counts:
            total_single_mismatched_counts += counts
        for (mismatch_pair, counts) in sorted_terminal_mismatch_counts:
            single_mismatched_ends_outfile.write(mismatch_pair[0] + '\t' + mismatch_pair[1] + '\t' + str(counts) + '\t' + str( float(counts)/total_single_mismatched_counts) + '\t' + str(total_single_mismatched_counts) + '\n')
        single_mismatched_ends_outfile.close()  

        three_prime_mismatch_seq_top_examples_outfile = open(filename_prefix + '_three_prime_mismatch_seq_top_examples.txt', 'w')
        num_examples = 100000
        three_prime_mismatch_seq_top_examples_outfile.write('three_prime_non_matching_seq\ttotal_three_prime_non_matching_counts\ttotal_counts_for_seq\tfeature_counts\tchrom\tcoord\tstrand\n')
        total_three_prime_non_matching_counts = total_soft_clipped_counts  + total_single_mismatched_counts 
        total_counts_for_seq = {}
        for three_prime_non_matching_seq in three_prime_non_matching_seq_to_locations_to_counts:
            total_counts_for_seq[three_prime_non_matching_seq] = 0
            for pA_coords in three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq]:
                pA_coords_counts = three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq][pA_coords]
                total_counts_for_seq[three_prime_non_matching_seq] += pA_coords_counts
        for three_prime_non_matching_seq in three_prime_non_matching_seq_to_locations_to_counts:
            sorted_pA_coord_counts = sorted(three_prime_non_matching_seq_to_locations_to_counts[three_prime_non_matching_seq].items(), key=operator.itemgetter(1), reverse = True)
            for idx in range(min(num_examples, len(sorted_pA_coord_counts) ) ):
                pA_coords, counts = sorted_pA_coord_counts[idx]
                chromosome, coordinate, strand = pA_coords
                soft_clipping_info = [three_prime_non_matching_seq, total_three_prime_non_matching_counts, total_counts_for_seq[three_prime_non_matching_seq], counts ]
                output = [str(e) for e in soft_clipping_info] + [str(e) for e in pA_coords]
                three_prime_mismatch_seq_top_examples_outfile.write('\t'.join(output) + '\n')
        three_prime_mismatch_seq_top_examples_outfile.close()    
            
