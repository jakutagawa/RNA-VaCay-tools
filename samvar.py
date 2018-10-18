#!/usr/bin/env python2
# Name: Jon Akutagawa
# Date: 18-08-14

'''
Adds simple variants into a BAM/SAM file. Reads in a text file with variant
locations, base change, allele frequency and mutation type. Outputs new BAM/SAM
file with variants.

Sample commands -
user$ python samvar.py -vf /scratch/jakutagawa/icgc/bams/test/test_snv.txt -ib /scratch/jakutagawa/icgc/bams/test/testregion_chr22.bam -is /scratch/jakutagawa/icgc/bams/test/testregion_chr22.sam -os /scratch/jakutagawa/icgc/bams/synthetic/testregion_chr22.with_variants.sam
user$ python samvar.py -vf /private/groups/brookslab/jakutagawa/variant_calling/synthetic_mutation_lists/October_2016_whitelist_2583.snv_mnv_indel.random_snp_only.txt -ib /scratch/jakutagawa/icgc/bams/normal/DO46933/PCAWG.764a33dc-dd34-11e4-8a0c-117cc254ba06.STAR.v1.bam -ob /scratch/jakutagawa/icgc/bams/synthetic/PCAWG.764a33dc-dd34-11e4-8a0c-117cc254ba06.STAR.v1.bam.with_variants.bam
'''
import sys
import random
import pysam
import re

from Bio import pairwise2

class Samvar :
    def __init__ (self, arguments):
        '''
        Initialize filename passed in and values needed for conversion.
        '''
        # create names for filenames
        self.variant_file = arguments.variant_file
        self.sam_input = arguments.sam_input
        self.bam_input = arguments.bam_input
        self.sam_output = arguments.sam_output
        self.bam_output = arguments.bam_output

        if self.sam_input:
            self.samfile = pysam.AlignmentFile(self.sam_input,"r")
            self.sam_header = self.samfile.header
        if self.bam_input:
            self.bamfile = pysam.AlignmentFile(self.bam_input,"rb")

        self.seed_num = arguments.seed

        # dictionary to match genome location to variants
        self.variants_dict = dict()

        # dictionary of read ids to variants
        self.reads_to_mutate = dict()

        # create dictionary matching location to list of reads
        self.reads_of_interest = dict()
        # dictionary matching read id/location to selected mutated read data
        self.mutated_reads = dict()


    def load_variants (self):
        '''
        Load variant file by calling read_variant_file and stores gene names and
        transcript ids into dictionary
        '''
        variant_counter = 0
        for variant in self.read_variant_file():
            variant_counter += 1
            #variant = (chromosome, pos1, pos2, vaf, base)
            chromosome = variant[0]
            pos1 = variant[1]
            pos2 = variant[2]
            location = chromosome + '_' + str(pos1) + '_' + str(pos2)
            self.variants_dict[location] = (variant)
            sys.stderr.write('Loaded '+ str(variant) + ' \n')
        sys.stderr.write('Loaded '+ str(variant_counter) + ' variants total \n')


    def read_variant_file (self):
        '''
        Read variant file line and splits line into separate variables.
        '''
        with open(self.variant_file) as fileH:
            for line in fileH:
                split_line = line.rstrip().split('\t')
                chromosome = split_line[0]
                pos1 = int(split_line[1])
                pos2 = int(split_line[2])
                vaf = float(split_line[3])
                try:
                    base = split_line[4]
                except IndexError:
                    base_dict = {1:'A',2:'T',3:'C',4:'G'}
                    random_base = base_dict[random.randint(1,4)]
                    base = random_base
                try:
                    type = split_line[5]
                except IndexError:
                    type = 'SNP'

                yield (chromosome,pos1,pos2,vaf,base,type)

    def count_and_shuffle_reads_with_variants (self):
        '''
        Uses variants/VAF to shuffle reads that have been successfully edited
        and returns a dictionary of shuffled reads
        '''


        for location, variant in self.variants_dict.items():
            complete_read_list = list()
            chromosome = variant[0]
            var_pos1 = variant[1]
            var_pos2 = variant[2]
            base = variant[4]
            vaf = variant[3]
            region = self.bamfile.fetch(chromosome,var_pos1,var_pos2+1)
            read_count = 0
            for read in region:
                ref_positions = read.get_reference_positions(full_length = True)
                if var_pos1 in ref_positions:
                    read_count += 1
                    complete_read_list.append((read.query_name, read.reference_start))

            reads_desired = int(vaf * read_count)
            random.seed(self.seed_num)
            shuffled_ids = list(complete_read_list)
            random.shuffle(shuffled_ids)
            shuffled_ids = shuffled_ids[:reads_desired]

            for read in shuffled_ids:
                try:
                    self.reads_to_mutate[read].append((var_pos1,var_pos2,base))
                except KeyError:
                    self.reads_to_mutate[read] = [(var_pos1,var_pos2,base)]

        sys.stderr.write('Randomly selected reads for mutation \n')

    def mutate_read (self, read):
        '''
        Add mutations to read. Calls edit_sequence and edit_cigar to make
        necessary changes.
        '''
        read_id = (read.query_name, read.reference_start)
        ref_positions = read.get_reference_positions(full_length = True)
        ref_seq = read.get_reference_sequence()
        old_seq = read.query_sequence
        old_query_seq = read.query_sequence
        ref_cigar = read.cigartuples
        old_cigar = ref_cigar
        old_md = read.get_tag('MD')
        new_md_tag = ''
        query_pos = (read.query_alignment_start, read.query_alignment_end)
        #new_clip_lengths = [0,0]

        #if len(self.reads_to_mutate[read_id]) > 2:
        #    print (read_id)
        #print (str(read_id) + ' has ' + str(len(self.reads_to_mutate[read_id])) + ' mutations')
        #print (old_cigar)

        for mutation in self.reads_to_mutate[read_id]:
            #print (mutation)
            var_pos1 = mutation[0]
            var_pos2 = mutation[1]
            new_base = mutation[2]

            #new_seq = self.edit_sequence(old_seq,new_base,var_pos1,ref_positions)
            #old_seq = new_seq
            new_query_seq, new_cigar = self.edit_seq_and_cigar(ref_seq,old_query_seq,old_cigar,new_base,var_pos1,ref_positions)
            old_query_seq = new_query_seq
            old_cigar = new_cigar
        #print(new_cigar)
        #new_cigar = new_attributes[1]

        #print (old_query_seq)
        #print (old_cigar)
        #print (old_md)

        #new_cigar = self.edit_cigar(old_cigar, read.query_alignment_sequence, query_pos, new_base, var_pos1, ref_positions)
        #old_cigar = new_cigar


        new_md_tag = self.edit_md_tag(old_md, ref_cigar, new_cigar, ref_seq, new_query_seq)
        new_nm_tag = self.edit_nm_tag(new_md_tag)

        #ref_seq = read.get_reference_sequence()
        #print (read.cigarstring)
        #print (read.get_reference_sequence())
        #print (read.query_sequence)
        #print (new_seq)
        #print (old_seq)
        #print (new_attributes[0])
        #new_seq = new_attributes[0]
        #new_cigar = new_attributes[1]
        #new_md = new_attributes[2]
        #if read.cigarstring != '101M':
        #    print (read.get_tag('MD'))
        #    print (read.get_aligned_pairs(with_seq=True))
        """
        if read.query_name == 'HISEQ4_0118:2:1316:1852:9553#GGCTAC':
            print (str(read.query_name) + ' added ' + str(len(self.reads_to_mutate[read_id])) + ' mutations')
            print (ref_seq)
            print (read.query_sequence)
            print (read.cigarstring)
            print (read.cigartuples)
            print (new_seq)
            '''
            print (new_cigar)
            gen_cigar = self.generate_cigar (ref_seq, new_seq, True)
            if ("I" in gen_cigar) or ("D" in gen_cigar):
                print (self.create_pairwise_alignment(ref_seq, new_seq, True))
                print (gen_cigar)
            else:
                print (gen_cigar)
            '''
        """
        qualities_copy = read.query_qualities
        read.query_sequence = new_query_seq
        read.query_qualities = qualities_copy
        read.cigartuples = new_cigar
        read.set_tag('MD', new_md_tag, value_type='Z')
        read.set_tag('NM', new_nm_tag, value_type='i')
        return read


    def edit_sequence (self, seq, new_base, var_pos, ref_positions):
        '''
        Add new base(s) to sequence
        '''
        if var_pos in ref_positions:
            var_index = ref_positions.index(var_pos)
            return (seq[:var_index] + new_base + seq[var_index+1:])
        else:
            return seq

    def create_pairwise_alignment(self, ref_seq=None, query_seq=None, local=True):
        '''
        Create alignment table from two strings using pairwise2 package
        :param ref_seq: string for reference squence
        :param query_seq: string for query sequence
        :param local: if true do local alignment else do global
        :return: dictionary with 'reference' and 'query' accessible alignment sequences
        '''
        assert ref_seq is not None, "Must set reference sequence"
        assert query_seq is not None, "Must set query sequence"

        if local:
            alignments = pairwise2.align.localms(ref_seq.upper(), query_seq.upper(), 1, -10, -10, -10,
                                                 one_alignment_only=True)

        else:
            alignments = pairwise2.align.globalms(ref_seq.upper(), query_seq.upper(), 2, -0.5, -1, -0.3,
                                                  one_alignment_only=True)
        # print(format_alignment(*alignments[0]))
        return {'reference': alignments[0][0], 'query': alignments[0][1]}

    def generate_cigar (self, ref_seq=None, query_seq=None, local=True):
        '''
        Generate CIGAR from alignment of two reads
        :param ref_seq: string for reference squence
        :param query_seq: string for query sequence
        :param local: bool for local or global alignment
        '''
        assert ref_seq is not None, "Must set reference sequence"
        assert query_seq is not None, "Must set query sequence"

        alignment = self.create_pairwise_alignment(ref_seq=ref_seq, query_seq=query_seq, local=local)
        final_str = str()
        # CIGAR = {"M": 0, "I": 0, "D": 0, "N": 0, "S": 0, "H": 0, "P": 0, "=": 0, "X": 0}
        current_op = str()
        op_count = 0
        for ref_seq, query_seq in zip(alignment['reference'], alignment["query"]):
            if ref_seq == query_seq:
                # matches
                if current_op == "M":
                    op_count += 1
                    current_op = "M"
                else:
                    final_str += str(op_count) + current_op
                    current_op = "M"
                    op_count = 1
            elif ref_seq == '-':
                # soft clipped sequences
                if current_op == "S":
                    op_count += 1
                    current_op = "S"
                elif current_op == str():
                    final_str += str(op_count) + current_op
                    current_op = "S"
                    op_count = 1
                # insertions
                elif current_op == "I":
                    op_count += 1
                    current_op = "I"
                else:
                    final_str += str(op_count) + current_op
                    current_op = "I"
                    op_count = 1
            elif query_seq == '-':
                # deletions
                if current_op == "D":
                    op_count += 1
                    current_op = "D"
                else:
                    final_str += str(op_count) + current_op
                    current_op = "D"
                    op_count = 1
            else:
                # mismatches
                if current_op == "S":
                    op_count += 1
                    current_op = "S"
                else:
                    final_str += str(op_count) + current_op
                    current_op = "S"
                    op_count = 1
        if current_op == "I":
            final_str += str(op_count) + "S"
        else:
            final_str += str(op_count) + current_op
        # remove initial zero
        return final_str[1:]

    def edit_nm_tag (self, md_tag):
        '''
        Count the number of mismatches in a MD tag to calculate the NM tag
        '''
        nm_count = sum(md_piece.isalpha() for md_piece in md_tag)
        return nm_count

    def edit_md_tag (self, md_tag, ref_cigar, new_cigar, ref_seq, new_query_seq):
        '''
        Edit existing MD tag based on mutation location and return updated MD
        tag
        '''
        split_md = re.split("(\D+)", md_tag)
        bases = ['A','C','T','G']
        md_index = 0
        new_md = list()
        old_clip_lengths = [0,0]
        new_clip_lengths = [0,0]
        clip_dif = [0,0]
        clipped_ref_seq = ref_seq
        clipped_query_seq = new_query_seq

        # check if ref and new cigar has a clip
        if len(ref_cigar) > 1 and len(new_cigar) > 1:

            # check for leading clip
            if ref_cigar[0][0] == 4 and new_cigar[0][0] == 4:
                old_clip_lengths[0] = ref_cigar[0][1]
                new_clip_lengths[0] = new_cigar[0][1]
                clip_dif[0] = new_cigar[0][1] - ref_cigar[0][1]

            # check for ending clip
            if ref_cigar[-1][0] == 4 and new_cigar[-1][0] == 4:
                old_clip_lengths[1] = ref_cigar[-1][1]
                new_clip_lengths[1] = new_cigar[-1][1]
                clip_dif[1] = new_cigar[-1][1] - ref_cigar[-1][1]

            # truncate sequences if necessary
            if clip_dif[0] == 0 and clip_dif[1] == 0:
                # check if front clip changes
                if new_clip_lengths[0]:
                    clipped_query_seq = clipped_query_seq[new_clip_lengths[0]:]
                # check if end clip changes
                if new_clip_lengths[1]:
                    clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]

            elif clip_dif[0] > 0:
                if new_clip_lengths[0]:
                    clipped_ref_seq = ref_seq[new_clip_lengths[0]-ref_cigar[0][1]:]
                    clipped_query_seq = clipped_query_seq[new_clip_lengths[0]:]
                if new_clip_lengths[1] and clip_dif[1]:
                    clipped_ref_seq = ref_seq[:-new_clip_lengths[1]+ref_cigar[-1][1]]
                    clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]
                elif new_clip_lengths[1] and clip_dif[1] == 0:
                    clipped_ref_seq = ref_seq[:-new_clip_lengths[1]]
                    clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]

            elif clip_dif[1] > 0:
                if new_clip_lengths[1]:
                    clipped_ref_seq = ref_seq[:-new_clip_lengths[1]+ref_cigar[-1][1]]
                    clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]
                if new_clip_lengths[0] and clip_dif[0]:
                    clipped_ref_seq = ref_seq[new_clip_lengths[0]-ref_cigar[0][1]:]
                    clipped_query_seq = clipped_query_seq[new_clip_lengths[0]:]
                elif new_clip_lengths[0] and clip_dif[0] == 0:
                    clipped_ref_seq = ref_seq[:-new_clip_lengths[1]]
                    clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]

        # check if soft clip gained at front or end
        elif len(ref_cigar) == 1 and len(new_cigar) > 1:
            # leading soft clip
            if new_cigar[0][0] == 4:
                new_clip_lengths[0] = new_cigar[0][1]
                clip_dif[0] = new_cigar[0][1]

            # check for ending clip
            if new_cigar[-1][0] == 4:
                new_clip_lengths[1] = new_cigar[-1][1]
                clip_dif[1] = new_cigar[-1][1]

            # clips sequences based on new soft clips
            if new_cigar[0][0] == 4:
                clipped_ref_seq = ref_seq[new_clip_lengths[0]:]
                clipped_query_seq = clipped_query_seq[new_clip_lengths[0]:]
            if new_cigar[-1][0] == 4:
                clipped_ref_seq = ref_seq[:-new_clip_lengths[1]]
                clipped_query_seq = clipped_query_seq[:-new_clip_lengths[1]]

        # find mismatching bases between sequences
        if clipped_ref_seq == clipped_query_seq:
            new_md = str(len(clipped_ref_seq))

        else:
            new_md = ''
            base_count = 0
            for index, base in enumerate(clipped_ref_seq):
                if base == clipped_query_seq[index]:
                    base_count += 1
                else:
                    new_md += str(base_count)
                    new_md += base
                    base_count = 0
            if base_count != 0:
                new_md += str(base_count)

        return new_md


    def edit_seq_and_cigar (self, ref_seq, query_seq, cigartuples, new_base, var_pos, ref_positions):
        '''
        Edit string, cigar string (mostly soft clips), and MD tag if changes
        are necessary
        '''
        if var_pos in ref_positions:
            #print (query_seq)
            var_index = ref_positions.index(var_pos)
            new_seq = (query_seq[:var_index] + new_base + query_seq[var_index+1:])
            new_cigar = ''
            clip_lengths = [0,0]

            # check for soft clips
            if len(cigartuples) > 1:
                # leading soft clip
                if cigartuples[0][0] == 4:
                    clip_lengths[0] = cigartuples[0][1]
                # ending soft clip
                if cigartuples[-1][0] == 4:
                    clip_lengths[1] = cigartuples[-1][1]


            # check if cigar is a perfect match - '101M'
            if cigartuples == [(0,101)]:
                # edge cases
                if var_index == 0:
                    new_cigar = [(4,1),(0,100)]
                elif var_index == 100:
                    new_cigar = [(0,100),(4,1)]
                else:
                    new_cigar = cigartuples

                return (new_seq, new_cigar)

            else:
                cigar_index = 0
                # loop through each cigar op tuple
                for op_index, op_tuple in enumerate(cigartuples):
                    cigar_index += op_tuple[1]
                    # convert M to S at the beginning
                    if op_index == 0 and var_index == cigar_index:
                        if op_tuple[0] == 4 and cigartuples[1][0] == 0:
                            first_tuple = list(cigartuples[0])
                            second_tuple = list(cigartuples[1])
                            first_tuple[1] += 1
                            second_tuple[1] -= 1
                            new_cigar = [tuple(first_tuple),tuple(second_tuple)] + cigartuples[2:]
                            return (new_seq, new_cigar)

                    # convert M to S at the end
                    elif op_index == len(cigartuples) - 2 and var_index == cigar_index:
                        if op_tuple[0] == 0 and cigartuples[-1][0] == 4:
                            last_tuple = list(cigartuples[-1])
                            penultimate_tuple = list(cigartuples[-2])
                            last_tuple[1] += 1
                            penultimate_tuple[1] -= 1
                            new_cigar = cigartuples[:-2] + [tuple(penultimate_tuple),tuple(last_tuple)]
                            return (new_seq, new_cigar)

                    elif var_index < cigar_index:
                        return (new_seq, cigartuples)

                return new_seq, cigartuples


            return (new_seq,cigartuples)
        else:
            return (query_seq,cigartuples)

    def clean_cigar (self, cigar):
        '''
        Remove 0S and 0M sections from cigar tuples in mutated reads and merge
        sections if necessary
        '''
        new_cigar = ''
        # grab index position of 0S or 0M
        first_element = [i[1] for i in cigar]

        if 0 in first_element:
            zero_index = first_element.index(0)
            cigar_pieces = len(cigar)
            new_cigar = ''

            cigar[zero_index]
            if zero_index-1 >= 0 and zero_index + 1 <= len(cigar):
                if cigar[zero_index+1][0] == cigar[zero_index-1][0]:
                    merged_element = (cigar[zero_index+1][0],cigar[zero_index+1][1]+cigar[zero_index-1][1])
                    new_cigar = [merged_element] + cigar[zero_index+2:]
                    return new_cigar

            else:
                sys.stderr.write ('mismatching cigar elements \n')
                return cigar

        else:
            return cigar

    def output_reads_to_bam (self):
        '''
        Loops through reads and outputs reads to new bam, mutating reads as
        necessary
        '''
        total_read_count = 0
        mutated_read_count = 0
        with pysam.AlignmentFile(self.bam_output, "wb", header = self.bamfile.header) as outf:
            for read in self.bamfile.fetch():
                total_read_count += 1
                if (total_read_count % 10000) == 0:
                    sys.stderr.write(str(total_read_count) + ' reads written \n')
                read_id = (read.query_name, read.reference_start)
                if read_id not in self.reads_to_mutate.keys():
                    outf.write(read)
                else:
                    mutated_read_count += 1
                    mutated_read = self.mutate_read(read)
                    outf.write(mutated_read)
        sys.stderr.write('New bam file: ' + str(self.bam_output) + '\n')
        sys.stderr.write('Total reads: ' + str(total_read_count) + '\n')
        sys.stderr.write('Mutated reads: ' + str(mutated_read_count) + '\n')

    def sort_and_index_bam (self):
        '''
        Sort and index output bam
        '''
        sorted_output = self.bam_output[:-4] + ".sorted.bam"
        pysam.sort("-o",sorted_output,self.bam_output)
        pysam.index(sorted_output)

        sys.stderr.write('New sorted bam file: ' + str(sorted_output) + '\n')


    def load_reads (self):
        '''
        Deprecated. Fetch BAM region associated with each variant and saves to a dictionary
        matching location to reads
        '''
        for location, variant in self.variants_dict.items():
            chromosome = variant[0]
            var_pos1 = variant[1]
            var_pos2 = variant[2]
            new_base = variant[4]

            region = self.bamfile.fetch(chromosome,var_pos1,var_pos2+1)
            read_count = self.bamfile.count(chromosome,var_pos1,var_pos2+1)
            #for read in region:
            #    read_count += 1
            print(read_count)
            random_line = list(range(read_count))
            print (random_line)
            random.shuffle(random_line)
            print (random_line)

            read_count = 0
            read_count2 = 0
            for x in region:

                #if (x.cigarstring != '101M'):
                var_index = var_pos1 - x.reference_start
                #print str(x.get_reference_name())
                '''
                print (x.cigarstring)
                print (x.cigar)
                print (x.reference_name)
                print (x.reference_start)


                print str(x.get_blocks())
                print (x.query_name)
                '''

                ref_positions = x.get_reference_positions(full_length = True)
                #clipped_ref_positions = x.get_reference_positions()
                #print str(x.get_reference_positions())

                if var_pos1 in ref_positions:
                    read_count2+=1


                new_seq = self.edit_sequence(x.query_sequence,new_base,var_pos1,ref_positions)
                '''
                if x.query_name == 'HISEQ4_0118:2:2116:2109:56685#GGCTAC':
                    print str(x.get_reference_positions(full_length = True))
                    print (x.cigarstring)
                    print (x.query_alignment_sequence)
                    print (x.query_sequence)
                    print (new_seq)
                '''

                #print (x.query_alignment_start)
                #print (x.query_alignment_end)
                query_pos = (x.query_alignment_start, x.query_alignment_end)


                #print (x)

                #print (x.query_qualities)

                if new_seq != x.query_sequence:
                    read_count += 1
                    qualities_copy = x.query_qualities
                    x.query_sequence = new_seq
                    x.query_qualities = qualities_copy
                    new_cigar = self.edit_cigar(x.cigartuples, x.query_alignment_sequence, query_pos, new_base, var_pos1, ref_positions)
                    x.cigartuples = new_cigar
                    #if x.query_name == 'HISEQ4_0118:2:2116:2109:56685#GGCTAC':
                    #    print str(x.get_reference_positions(full_length = True))


                    #print (x.cigarstring)
                    #self.all_reads_of_interest.append()
                    try:
                        self.reads_of_interest[location].append([x.query_name,x.cigarstring,x.reference_start,x.query_sequence])
                        #self.reads_of_interest[location].append(x)
                    except KeyError:
                        self.reads_of_interest[location] = [[x.query_name,x.cigarstring,x.reference_start,x.query_sequence]]
                        #self.reads_of_interest[location] = [x]
            print (read_count)
            print (read_count2)
        print str((self.reads_of_interest).keys())

    def print_bam_file (self):
        '''
        Deprecated. Print info from a pileup.
        '''
        iter = self.bamfile.pileup('22', 38062455, 38062456)
        column_count = 0
        for pileupcolumn in iter:
            #print (str(x.get_query_names()))
            if pileupcolumn.n > 80:
                column_count += 1
                print ("\ncoverage at base %s = %s" %
                    (pileupcolumn.pos, pileupcolumn.n))
                read_count = 0
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        read_count += 1
                        print ('\t %s base in read %s = %s at pos %s' %
                            (read_count, pileupread.alignment.query_name,
                            pileupread.alignment.query_sequence[pileupread.query_position],pileupread.query_position))
        sys.stderr.write('column count: ' + str(column_count) + '\n')

    def mutate_reads (self):
        '''
        Deprecated. Mutate read by inserting new base(s) to sequence and modifying cigar
        string as necessary
        '''
        for location, read_list in self.reads_of_interest.items():
            for read in read_list:

                split_location = location.split('_')
                var_pos1 = int(split_location[1])
                var_pos2 = int(split_location[2])
                new_base = self.variants_dict[location][4]

                '''
                id = read[0]
                old_cigar = read[1]
                ref_start = int(read[2])
                '''
                id = read.query_name
                old_cigar = read.cigar
                ref_start = int(read.reference_start)

                print ('old cigar is ' + old_cigar)
                #print (read)

                var_index = var_pos1 - ref_start
                new_read = list()
                # check if cigar is a perfect match
                if old_cigar == '101M':
                    #pass


                    seq = read[3]
                    #print (read)
                    #print (var_index)
                    #split_cigar = self.split_cigar(cigar)
                    #print (split_cigar)
                    #print (seq[var_index])
                    #print (new_base)

                    if var_index == 0:
                        new_cigar = str(var_index + 1) + 'S' + str(101 - var_index - 1) + 'M'
                    elif var_index == 100:
                        new_cigar = str(var_index) + 'M' + str(101 - var_index) + 'S'
                    else:
                        new_cigar = str(var_index) + 'M1S' + str(101 - var_index - 1) + 'M'
                    #print (new_cigar)
                    #print (seq)
                    new_seq = self.edit_sequence(seq,new_base,var_index)
                    #print (read)
                    new_read = list(read)
                    new_read[1] = new_cigar
                    new_read[3] = new_seq
                    print (read)
                    print (new_read)
                else:

                    print (var_index)
                    split_cigar = self.split_cigar(old_cigar)
                    cigar_pieces = len(split_cigar)/2
                    cigar_sum = 0
                    cigar_sum_no_n = 0
                    #for val in range(0,cigar_pieces):


                    seq_counter = 0
                    new_cigar = ''
                    for val in range(0,cigar_pieces):
                        cigar_sum += split_cigar[2*val]
                        if split_cigar[2*val+1] != 'N':
                            cigar_sum_no_n += int(split_cigar[2*val])
                        old_counter = seq_counter
                        seq_counter += split_cigar[2*val]
                        if var_index < seq_counter and var_index >= old_counter and split_cigar[2*val+1] != 'N':
                            print (var_index)
                            #if var_index == 0:
                            #    new_cigar += str(var_index + 1)
                            #elif var_index == cigar_sum - 1:
                            #    new_cigar += str(var_index)
                            if var_index > 0:
                                new_cigar += str(var_index - old_counter) + split_cigar[2*val+1]
                            #elif var_index + 1 == seq_counter:
                            #    new_cigar += str(var_index - old_counter) +

                            if split_cigar[2*val+1] == 'M':
                                new_cigar += '1S'
                            elif split_cigar[2*val+1] == 'S':
                                #if seq[var_index] ==
                                new_cigar += '1M'

                            new_cigar += str(seq_counter - var_index - 1) + split_cigar[2*val+1]
                        else:
                            new_cigar += str(split_cigar[2*val])
                            new_cigar += split_cigar[2*val+1]
                    #new_cigar = str(var_index)
                    print ('new cigar is ' + new_cigar)
                    #new_split_cigar = self.split_cigar(new_cigar)
                    merged_new_cigar = self.merge_cigar(new_cigar)
                    print (merged_new_cigar)

                    #if 'N' in split_cigar:
                    #    var_index = var_pos1 - ref_start
                    #else:
                    #print (split_cigar)

                    #print (split_cigar)
                    #print (cigar_pieces)
                    '''
                    for val in range(0,cigar_pieces):
                        print(split_cigar[2*val])
                        if split_cigar[2*val+1] == 'S':
                            elif split_cigar[2*val+1] == ''


                    else:
                        for val in range(0,cigar_pieces):
                            print(split_cigar[2*val])
                            if split_cigar[2*val+1] == 'M' and split_cigar[2*val] > var_index:
                                elif split_cigar[2*val+1] == ''
                        if var_index < split_cigar[0]:
                            new_cigar = str(var_index) + 'M' + str(var_index + 1) + 'S' + str(101 - var_index - 1)

                    '''
                    #seq = read[3]
                    #var_index = var_pos1 - ref_start
                    #print(seq[:var_index])
                    '''
                    ref_start = int(read[2])
                    var_index = var_pos1 - ref_start
                    seq = read[3]
                    print (read)
                    print (var_index)
                    print (seq)
                    print (self.split_cigar(cigar))
                    print (seq[var_index])
                    '''
                if new_read:
                    self.mutated_reads[(id,ref_start)] = new_read
        print (self.mutated_reads)

    def merge_cigar (self, cigar):
        '''
        Deprecated. Remove 0S sections from cigar strings in mutated reads
        '''
        split_cigar = self.split_cigar(cigar)
        new_cigar = ''
        if 0 in split_cigar:
            cigar_pieces = len(split_cigar)/2
            new_cigar = ''
            for val in range(0,cigar_pieces):
                if split_cigar[2*val] == 0:
                    # check if 0S at start
                    if val == 0:
                        new_cigar = split_cigar[2:]
                    # check if 0S at end
                    elif val == cigar_pieces - 1:
                        new_cigar = split_cigar[:-2]
                    else:
                        # check if operations match
                        if split_cigar[2*val-1] == split_cigar[2*val+3]:

                            new_cigar = split_cigar[:2*val-2]
                            #print ('new_cigar is '+ str(new_cigar))
                            new_cigar += [str(split_cigar[2*val-2]+split_cigar[2*val+2]), split_cigar[2*val-1]]
                            #print ('new_cigar is '+ str(new_cigar))
                            try:
                                new_cigar += split_cigar[2*val+4:]
                            except IndexError:
                                pass
                            #print ('new_cigar is '+ str(new_cigar))
                        else:
                            new_cigar = split_cigar[:2*val] + split_cigar[2*val+2:]


                #print ('new_cigar is '+ str(new_cigar))
        else:
            new_cigar = split_cigar
        #print (split_cigar)

        merged_cigar = ''.join([str(i) for i in new_cigar])
        #print (merged_cigar)
        return merged_cigar

    def split_cigar (self, cigar):
        '''
        Deprecated. Split cigar string into a list of integers and characters
        '''
        split_cigar = filter (None, re.split(r'(\d+)',cigar))

        for index, piece in enumerate(split_cigar):
            if index % 2 == 0:
                split_cigar[index] = int(piece)

        return split_cigar

    def convert_to_sam (self):
        '''
        Incomplete. Convert input BAM to SAM format
        '''
        pysam.view()

    def output_reads_to_sam (self):
        '''
        Deprecated and incomplete. Read through and write a sam file.
        '''
        for read in self.read_sam_file():
            mut_read_loc = (read[0],read[2])
            if mut_read_loc in self.mutated_reads.keys():
                mutated_read = self.mutated_reads[mut_read_loc]
                sys.stdout.write(mutated_read[0] + '\t' + mutated_read[1])
            else:
                sys.stdout.write()

    def read_sam_file (self):
        '''
        Deprecated. Read sam file line by line and splits line into separate variables.
        '''
        with open(self.sam_input) as fileH:
            line = fileH.readline()
            while line.startswith('@'):
                line = fileH.readline()
                sys.stderr.write(line + '\n')
            for line in fileH:
                split_line = line.rstrip().split('\t')
                #sys.stderr.write(str(split_line) + '\n')
                qname = split_line[0]
                flag = split_line[1]
                chr = split_line[2]
                pos = split_line[3]
                mapq = split_line[4]
                cigar = split_line[5]
                mrnm = split_line[6]
                mate_pos = split_line[7]
                tlen = split_line[8]
                sequence = split_line[9]
                qual = split_line[10]
                tags = split_line[11]

                yield (qname,flag,chr,pos1,pos2,cigar,sequence)




class CommandLine() :
    '''
    modified from David Bernick

    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    my_command_line.args is a dictionary which includes each of the available command line arguments as
    my_command_line.args['option']
    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - Specify the <input >output files and conditions',
                                             epilog = 'Program epilog - parameters of infile and outfile must be given.',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s -vf variant_file -is sam_input -os sam_output'
                                             )
        self.parser.add_argument('-vf','--variantFile', dest='variant_file',
                                 action='store', default = '', type=str,
                                 required=True,
                                 help='get variant list from file')
        self.parser.add_argument('-ib','--inputBAM', dest='bam_input',
                                 action='store', default = '', type=str,
                                 required=True,
                                 help='input BAM file location')
        self.parser.add_argument('-is','--inputSAM', dest='sam_input',
                                 action='store', default = '', type=str,
                                 required=False,
                                 help='input SAM file location')
        self.parser.add_argument('-os','--outputSAM', dest='sam_output',
                                 action='store', default = '', type=str,
                                 required=False,
                                 help='output SAM file location')
        self.parser.add_argument('-ob','--outputBAM', dest='bam_output',
                                 action='store', default = '', type=str,
                                 required=True,
                                 help='output BAM file location')
        self.parser.add_argument('-rs','--randomSeed', dest='seed',
                                 action='store', type=int, default=10,
                                 help='specify seed for randomizing reads')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

def main(my_command_line=None):
    '''
    Instantiate Samvar class. Get mutation BED and gene names from sys
    input filenames. Can output to a sys output file.
    '''
    my_command_line = CommandLine(my_command_line)
    my_samvar = Samvar(my_command_line.args)

    my_samvar.load_variants()

    #my_samvar.read_sam_file()
    #my_samvar.load_reads()
    my_samvar.count_and_shuffle_reads_with_variants()
    my_samvar.output_reads_to_bam()
    #my_samvar.sort_and_index_bam()
    #print(my_samvar.generate_cigar('GGTACCAGTTTAGGTTCCTAAGTAATAGTGACCCTTTCACGTCCTGGAGCCCGAGTGGACCAATCGGAAGCCTAAGTGACGATGACCCTCGCATGCCCTAG','GGTACAAGTTTAGGTTCCTAAGTAATAGTGACCCTTTCACGTCCTGGAGCCCGAGTGGACCAATCGGAAGCCCAAGTGACGATGACCCTCGCATGCCCTAG',True))


if __name__ == "__main__":
    main()
