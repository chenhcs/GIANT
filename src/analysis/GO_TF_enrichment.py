#!/usr/bin/env python
#
# This GIANT code is adapted from:
# Copyright (c) 2010-2018 Haibao Tang
# Licensed under the under the BSD 2-Clause "Simplified" License

import sys
import os.path as op
#from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc
from goatools.godag.consts import RELATIONSHIP_LIST
from goatools.multiple_testing import Methods
from goatools.pvalcalc import FisherFactory
import os
import glob
import argparse
import numpy as np
from pyensembl import EnsemblRelease

from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag

sys.path.insert(0, op.join(op.dirname(__file__), "."))

class GoeaCliArgs:
    """Extracts arguments from the command-line."""

    def __init__(self):
        self.args = self._init_args()

    def _init_args(self):
        """Get enrichment arg parser."""

        #pylint: disable=invalid-name
        p = argparse.ArgumentParser(__doc__,
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        p.add_argument('filenames', type=str, nargs=3,
                       help='data/study data/population data/association')
        p.add_argument('--annofmt', default=None, type=str,
                       help=('Annotation file format. '
                             'Not needed if type can be determined using filename'),
                       choices=['gene2go', 'gaf', 'gpad', 'id2gos'])
        p.add_argument('--taxid', default=9606, type=int,
                       help="When using NCBI's gene2go annotation file, specify desired taxid")
        p.add_argument('--alpha', default=0.05, type=float,
                       help='Test-wise alpha for multiple testing')
        p.add_argument('--pval', default=.05, type=float,
                       help=('Only print results with uncorrected p-value < PVAL. '
                             'Print all results, significant and otherwise, '
                             'by setting --pval=1.0'))
        p.add_argument('--pval_field', type=str,
                       help='Only print results when PVAL_FIELD < PVAL.')
        p.add_argument('--outfile', default=None, type=str,
                       help='Write enrichment results into xlsx or tsv file')
        p.add_argument('--ns', default='BP,MF,CC', type=str,
                       help='Limit GOEA to specified branch categories. '
                            'BP=Biological Process; '
                            'MF=Molecular Function; '
                            'CC=Cellular Component')
        p.add_argument('--id2sym', default=None, type=str,
                       help='ASCII file containing one geneid and its symbol per line')
        p.add_argument('--sections', default=None, type=str,
                       help=('Use sections file for printing grouped GOEA results. '
                             'Example SECTIONS values:\n'
                             'goatools.test_data.sections.gjoneska_pfenning \n'
                             'goatools/test_data/sections/gjoneska_pfenning.py \n'
                             'data/gjoneska_pfenning/sections_in.txt\n'))
        p.add_argument('--outfile_detail', type=str,
                       help=('Write enrichment results into a text file \n'
                             'containing the following information: \n'
                             '1) GOEA GO terms, grouped into sections \n\n'
                             '2) List of genes and ASCII art showing section membership \n'
                             '3) Detailed list of each gene and GO terms w/their P-values \n'))
        p.add_argument('--compare', dest='compare', default=False,
                       action='store_true',
                       help="the population file as a comparison group. if this "
                       "flag is specified, the population is used as the study "
                       "plus the `population/comparison`")
        p.add_argument('--ratio', dest='ratio', type=float, default=None,
                       help="only show values where the difference between study "
                       "and population ratios is greater than this. useful for "
                       "excluding GO categories with small differences, but "
                       "containing large numbers of genes. should be a value "
                       "between 1 and 2. ")
        p.add_argument('--prt_study_gos_only', default=False, action='store_true',
                       help=('Print GO terms only if they are associated with study genes. '
                             'This is useful if printng all GO results '
                             'regardless of their significance (--pval=1.0). '))
        p.add_argument('--indent', dest='indent', default=False,
                       action='store_true', help="indent GO terms")
        p.add_argument('--obo', default="go-basic.obo", type=str,
                       help="Specifies location and name of the obo file")
        p.add_argument('--no_propagate_counts', default=False, action='store_true',
                       help="Do not propagate counts to parent terms")
        # no -r:   args.relationship == False
        # -r seen: args.relationship == True
        p.add_argument('-r', '--relationship', action='store_true',
                       help='Propagate counts up all relationships,')
        # NO --relationships                -> None
        # --relationships part_of regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of           -> relationships=['part_of']
        # --relationships=part_of,regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of regulates -> NOT VALID
        p.add_argument('--relationships', nargs='*',
                       help=('Propagate counts up user-specified relationships, which include: '
                             '{RELS}').format(RELS=' '.join(RELATIONSHIP_LIST)))
        p.add_argument('--method', default="bonferroni,sidak,holm,fdr_bh", type=str,
                       help=Methods().getmsg_valid_methods())
        p.add_argument('--pvalcalc', default="fisher_scipy_stats", choices=FisherFactory.options.keys(),
                       help=str(FisherFactory()))
        p.add_argument('--min_overlap', default=0.7, type=float,
                       help="Check that a minimum amount of study genes are in the population")
        p.add_argument('--goslim', default='goslim_generic.obo', type=str,
                       help="The GO slim file is used when grouping GO terms.")
        p.add_argument('--ev_inc', type=str,
                       help="Include specified evidence codes and groups separated by commas")
        p.add_argument('--ev_exc', type=str,
                       help="Exclude specified evidence codes and groups separated by commas")
        p.add_argument('--ev_help', dest='ev_help', action='store_false',
                       help="Print all Evidence codes, with descriptions")
        p.add_argument('--ev_help_short', dest='ev_help_short', action='store_false',
                       help="Print all Evidence codes")
        # remove_goids: TBD
        #   None (Default) Remove a small number (~14) of broad GO IDs from the association
        #   True           Remove a slightly larger number of broad GO IDs (~100)
        #   False          Do not remove any broad GO IDs
        ## p.add_argument('--remove_goids', dest='remove_goids', default=None,
        ##                help="User-specified list of broad GO IDs to remove")

        if len(sys.argv) == 1:
            sys.exit(not p.print_help())
        self._prt_evidence_codes(set(sys.argv[1:]))
        args = p.parse_args()  # Namespace object from argparse
        self._adjust_relationships(args)
        self._check_input_files(args, p)
        return args

    @staticmethod
    def _prt_evidence_codes(args):
        if not {'--ev_help', '--ev_help_short'}.isdisjoint(args):
            print('\nEVIDENCE CODE HELP: --ev_exc --ev_inc')
            print('Use any of these group names, ')
            print('like Experimental or Similarity or Experimental,Similarity,')
            print('or evidence codes, like IEA or ISS,ISO,ISA in --ev_exc or --ev_inc:')
            obj = EvidenceCodes()
            if '--ev_help' in args:
                print('')
                obj.prt_details()
            if '--ev_help_short' in args:
                print('')
                obj.prt_summary_code()
            sys.exit(0)

    @staticmethod
    def _check_input_files(nspc, parser):
        """check filename args. otherwise if one of the 3 filenames is bad
        it's hard to tell which one"""
        if not len(nspc.filenames) == 3:
            parser.print_help()
            msg = """
      3 Expected files; Expected content: study population association",
      {} Actual   files: {}""".format(len(nspc.filenames), ' '.join(nspc.filenames))
            raise Exception(msg)
        for fin in nspc.filenames:
            if not os.path.exists(fin):
                return "*{}* does not exist".format(fin)
        return False

    @staticmethod
    def _adjust_relationships(args):
        """Adjust relationships for various user input"""
        # NO --relationships                -> None
        # --relationships part_of regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of           -> relationships=['part_of']
        # --relationships=part_of,regulates -> relationships=['part_of,regulates']
        # --relationships=part_of regulates -> NOT VALID
        if args.relationship:
            args.relationships = RELATIONSHIP_SET
        if args.relationships is not None:
            if len(args.relationships) == 1 and ',' in args.relationships[0]:
                args.relationships = args.relationships[0].split(',')
            args.relationships = set(args.relationships)
            chk_relationships(args.relationships)


def main():
    """Run gene enrichment analysis."""
    # Load study, population, associations. Run GOEA.
    obj = GoeaCliFnc(GoeaCliArgs().args)

    gene_lists = glob.glob(obj.args.filenames[0].split('/')[0] + '/cluster*')

    try:
        os.mkdir('/'.join(obj.args.filenames[0].split('/')[:-1]) + '/goea')
    except OSError as error:
        print(error)

    gene_lists.sort()
    for f in gene_lists:
        _study, _pop = obj.rd_files(f, obj.args.filenames[1])
        print(len(_pop), len(_study))
        obj.results_all = obj.objgoeans.run_study(_study)
        # Reduce results to significant results (pval<value)
        results_specified = obj.get_results()
        # Print results in a flat list
        outf = '/'.join(obj.args.filenames[0].split('/')[:-1]) + '/goea/'
        obj.args.outfile = outf + f.split('/')[-1]
        obj.prt_results(results_specified)

if __name__ == "__main__":
    main()
