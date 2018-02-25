#!/usr/bin/env python3
# Name: Jon Akutagawa
# Date: 18-02-23

"""
This script pulls random mutations from the PCAWG consensus variant call list to
generate a list of mutations in bed-like file for BAMSurgeon.

input: MAF file
output:
chr pos1    pos2    vaf base
22  234234  234237  0.25    A

python generate_random_mutations.py -t SNP -mf October_2016_whitelist_2583.snv_mnv_indel.maf.short >October_2016_whitelist_2583.snv_mnv_indel.maf.random
"""

import sys
import random

class RandomMutationListGenerator :
    def __init__ (self, arguments):
        """
        Initialize filename passed in and values needed for randomization.
        """
        # creates names for storing argument values
        self.maf_filename = arguments.maf_filename
        self.mut_type = arguments.mut_type
        self.seed_num = arguments.seed_num

        self.all_mutations = list()
        self.random_mutations = list()


    def load_maf (self):
        """
        Load all data into main list
        """
        counter = 0
        empty_counter = 0
        overlapping_counter = 0

        for chrom,pos1,pos2,mut_type,mutation,vaf in self.read_maf(self.maf_filename):
            counter += 1
            new_base = mutation.split('>')[1]
            mut_entry = (chrom,pos1,pos2,vaf,new_base)
            self.all_mutations.append(mut_entry)

        sys.stderr.write('loaded ' + str(counter) + ' mutations \n')

    def randomize (self):
        """
        Randomize list with user-defined seed
        """
        random.seed( self.seed_num )
        self.random_mutations = random.sample(self.all_mutations, 1000)

    def outputData (self):
        """
        Print data to sys output
        """
        for count, mutation in enumerate(self.random_mutations, 1):
            sys.stdout.write(str(mutation[0]) + '\t' + str(mutation[1]) + '\t')
            sys.stdout.write(str(mutation[2]) + '\t' + str(mutation[3]) + '\t')
            sys.stdout.write(str(mutation[4]) + '\n')

        sys.stderr.write('output ' + str(count) + ' random mutations \n')

    def read_maf (self, maf_filename):
        '''
        Read in MAF and parse line for necessary values.
        '''

        with open(maf_filename) as fileH:
            # read the header
            line = fileH.readline()
            header = line.rstrip()

            counter = 0
            for line in fileH:
                counter += 1
                split_line = line.rstrip().split('\t')
                chrom = split_line[1]
                pos1 = int(split_line[2])
                pos2 = int(split_line[3])
                mut_type = split_line[6]
                mutation = split_line[14]
                line_num = counter

                # some VAF values are split or missing
                vaf = split_line[28].split('|')[0]
                if vaf:
                    float_vaf = float(vaf)
                else:
                    vaf = 0

                # only yield if mutation type is correct and vaf not 0
                if mut_type in (self.mut_type) and vaf > 0:
                    yield chrom,pos1,pos2,mut_type,mutation,vaf




class CommandLine() :
    '''
    modified from David Bernick

    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']
    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - Specify the input files and conditions',
                                             epilog = 'Program epilog - parameters of infiles must be given.',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s -t SNP -mf mut.maf'
                                             )
        self.parser.add_argument('-t','--mutationType', dest='mut_type',
                                 action='store',type=str, required=True, nargs='+',
                                 help='mutation types desired')

        self.parser.add_argument('-mf','--mafFilename', dest='maf_filename',
                                 action='store',type=str, required=True,
                                 help='MAF filename')

        self.parser.add_argument('-s','--randomSeed', dest='seed_num',
                                 action='store',type=int, default=15,
                                 help='random seed number')
        """
        self.parser.add_argument('-of','--outputFile', dest='output_file',
                                 action='store',type=str, required=True,
                                 help='output data filename')
        """


        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(myCommandLine=None):
    """
    Instantiate RandomMutationListGenerator class. Load MAF. Output random mutations.
    """
    my_command_line = CommandLine(myCommandLine)
    random_mutation_list = RandomMutationListGenerator(my_command_line.args)

    random_mutation_list.load_maf()
    random_mutation_list.randomize()
    random_mutation_list.outputData()


if __name__ == "__main__":
    main()
