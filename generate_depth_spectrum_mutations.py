#!/usr/bin/env python3
# Name: Jon Akutagawa
# Date: 18-12-26

"""
This script pulls mutations from the PCAWG consensus variant call list to
generate a list of mutations in bed-like file for BAMSurgeon.

input: MAF file and BAM file
output:
chr pos1    pos2    vaf base    mutation_type
22  234234  234234  0.25    A   SNP

python generate_depth_spectrum_mutations.py -t SNP -mf /private/groups/brookslab/PCAWG/Oct2016_Freeze/October_2016_whitelist_2583.snv_mnv_indel.maf -n 300 -vc 5'UTR 3'UTR Missense_Mutation Nonsense_Mutation Nonstop_Mutation Silent Start_Codon_SNP >October_2016_whitelist_2583.snv_mnv_indel.random_snp_only.txt
"""

import sys
import random
import pysam

from .generate_random_mutations import RandomMutationListGenerator

class CommandLine() :
    """
    modified from David Bernick

    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']
    """

    def __init__(self, inOpts=None) :
        """
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        """
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - Specify the input files and conditions',
                                             epilog = 'Program epilog - parameters of infiles must be given.',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s -t SNP -mf mut.maf -n 200'
                                             )
        self.parser.add_argument('-t','--mutationType', dest='mut_type',
                                 action='store', type=str, required=True,
                                 nargs='+', help='mutation types desired')
        self.parser.add_argument('-vc','--variantClass', dest='variant_classes',
                                 action='store', type=str, required=False,
                                 nargs='+', help='variant classes desired')

        self.parser.add_argument('-mf','--mafFilename', dest='maf_filename',
                                 action='store', type=str, required=True,
                                 help='MAF filename')

        self.parser.add_argument('-rv','--randomVAF', dest='random_vaf',
                                 action='store_true', default=False,
                                 help='randomize VAFs')
        self.parser.add_argument('-s','--randomSeed', dest='seed_num',
                                 action='store',type=int, default=15,
                                 help='random seed number')
        self.parser.add_argument('-n','--mutationCount', dest='mut_num',
                                 action='store',type=int, default=1000,
                                 help='number of mutations desired')
        """
        self.parser.add_argument('-of','--outputFile', dest='output_file',
                                 action='store',type=str, required=True,
                                 help='output data filename')
        """


        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(my_command_line=None):
    """
    Instantiate RandomMutationListGenerator class. Load MAF. Select random
    mutations. Output random mutations to sys output.
    """
    my_command_line = CommandLine(my_command_line)
    random_mutation_list = RandomMutationListGenerator(my_command_line.args)

    random_mutation_list.load_maf()
    else:
        random_mutation_list.randomize()
    random_mutation_list.output_data()


if __name__ == "__main__":
    main()
