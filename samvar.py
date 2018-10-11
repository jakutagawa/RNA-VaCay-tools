#!/usr/bin/env python2
# Name: Jon Akutagawa
# Date: 18-08-14

"""
Adds simple variants into a BAM/SAM file. Reads in a text file with variant
locations, base change, allele frequency and mutation type. Outputs new BAM/SAM
file with variants.

Sample commands -
user$ python samvar.py -vf /scratch/jakutagawa/icgc/bams/test/test_snv.txt -ib /scratch/jakutagawa/icgc/bams/test/testregion_chr22.bam -is /scratch/jakutagawa/icgc/bams/test/testregion_chr22.sam -os /scratch/jakutagawa/icgc/bams/synthetic/testregion_chr22.with_variants.sam
user$ python samvar.py -vf /scratch/jakutagawa/icgc/bams/test/test_snv.txt -ib /scratch/jakutagawa/icgc/bams/test/testregion_chr22.bam -is /scratch/jakutagawa/icgc/bams/test/testregion_chr22.sam -os /scratch/jakutagawa/icgc/bams/synthetic/testregion_chr22.with_variants.sam
"""
import sys
import random
import pysam
import re

class Samvar :
    def __init__ (self, arguments):
        """
        Initialize filename passed in and values needed for conversion.
        """
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
        old_seq = read.query_sequence
        old_cigar = read.cigartuples
        query_pos = (read.query_alignment_start, read.query_alignment_end)

        for mutation in self.reads_to_mutate[read_id]:

            var_pos1 = mutation[0]
            var_pos2 = mutation[1]
            new_base = mutation[2]

            new_seq = self.edit_sequence(old_seq,new_base,var_pos1,ref_positions)
            old_seq = new_seq

            new_cigar = self.edit_cigar(old_cigar, read.query_alignment_sequence, query_pos, new_base, var_pos1, ref_positions)
            old_cigar = new_cigar


        qualities_copy = read.query_qualities
        read.query_sequence = new_seq
        read.query_qualities = qualities_copy
        read.cigartuples = new_cigar
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

    def edit_cigar (self, cigartuples, alignment_seq, query_pos, new_base, var_pos, ref_positions):
        '''
        Edit cigar sequence if changes necessary
        '''
        if var_pos in ref_positions:
            var_index = ref_positions.index(var_pos)
            new_cigar = ''

            # check if cigar is a perfect match - '101M'
            if cigartuples == [(0,101)]:
                # edge cases
                if var_index == 0:
                    new_cigar = [(4,1),(0,100)]
                elif var_index == 100:
                    new_cigar = [(0,100),(4,1)]
                else:
                    new_cigar = [(0,var_index),(4,1),(0,101 - var_index - 1)]
                return new_cigar
            else:
                if alignment_seq[var_index - query_pos[0]] == new_base:
                    new_op = 0
                else:
                    new_op = 4
                    #print(test)
                old_pos = 0
                for cigar_index, cigar_tuple in enumerate(cigartuples):
                    if var_index + 1 >= old_pos and var_index + 1 < old_pos + cigar_tuple[1]:
                        if var_index == 0:
                            first_cigar = list(cigartuples[0])
                            first_cigar[1] -= 1
                            #print ('first index')
                            new_cigar = [(new_op,1)] + [tuple(first_cigar)] + cigartuples[1:]

                        elif var_index == 100:
                            last_tuple = list(cigartuples[-1])
                            last_tuple[1] -= 1
                            #print ('last index')
                            new_cigar = cigartuples[:-1] + [tuple(last_tuple)] + [(new_op,1)]

                        else:
                            mid_tuple1 = list(cigartuples[cigar_index])
                            mid_tuple2 = list(cigartuples[cigar_index])
                            #print ('middle tuple')
                            mid_tuple1[1] = (var_index - old_pos)
                            mid_tuple2[1] -= (var_index - old_pos) + 1
                            new_cigar = cigartuples[:cigar_index] + [tuple(mid_tuple1)] + [(new_op,1)] + [tuple(mid_tuple2)] + cigartuples[cigar_index+1:]

                        cleaned_new_cigar = self.clean_cigar(new_cigar)
                        return cleaned_new_cigar
                    old_pos += cigar_tuple[1]

        else:
            return cigartuples

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
                """
                print (x.cigarstring)
                print (x.cigar)
                print (x.reference_name)
                print (x.reference_start)


                print str(x.get_blocks())
                print (x.query_name)
                """

                ref_positions = x.get_reference_positions(full_length = True)
                #clipped_ref_positions = x.get_reference_positions()
                #print str(x.get_reference_positions())

                if var_pos1 in ref_positions:
                    read_count2+=1


                new_seq = self.edit_sequence(x.query_sequence,new_base,var_pos1,ref_positions)
                """
                if x.query_name == 'HISEQ4_0118:2:2116:2109:56685#GGCTAC':
                    print str(x.get_reference_positions(full_length = True))
                    print (x.cigarstring)
                    print (x.query_alignment_sequence)
                    print (x.query_sequence)
                    print (new_seq)
                """

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

                """
                id = read[0]
                old_cigar = read[1]
                ref_start = int(read[2])
                """
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
                    """
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

                    """
                    #seq = read[3]
                    #var_index = var_pos1 - ref_start
                    #print(seq[:var_index])
                    """
                    ref_start = int(read[2])
                    var_index = var_pos1 - ref_start
                    seq = read[3]
                    print (read)
                    print (var_index)
                    print (seq)
                    print (self.split_cigar(cigar))
                    print (seq[var_index])
                    """
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
    """
    Instantiate Samvar class. Get mutation BED and gene names from sys
    input filenames. Can output to a sys output file.
    """
    my_command_line = CommandLine(my_command_line)
    my_samvar = Samvar(my_command_line.args)

    my_samvar.load_variants()

    #my_samvar.read_sam_file()
    #my_samvar.load_reads()
    my_samvar.count_and_shuffle_reads_with_variants()
    my_samvar.output_reads_to_bam()
    my_samvar.sort_and_index_bam()


if __name__ == "__main__":
    main()
