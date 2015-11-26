#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 14:14:33 2015

@author: sven
"""


import argparse
import textwrap
import sys
import CLQC
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='XiInfo',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        Scrip to create a distogram from cross-link identifications.
        Tested with HSA only so far.

        X-axes: Retention time
        Y-axes: summed intensity for MS1 scan

        ---------------------------------------------------------------------
        The program requires an input (mzml) and an output parameter (outfile)
        for the created image. The image can also bedisplayed directly.
        ---------------------------------------------------------------------
        '''))

    #mandatory parameters
    parser.add_argument('--in_links', metavar='in_links', type=str, nargs='+',
                        help='an integer for the accumulator')
    parser.add_argument('--in_fasta', metavar='in_fasta', type=str, nargs=1,
                        help='FASTA single proteine sequence file.')
    parser.add_argument('--out_path', metavar='out_path', type=str, nargs=1,
                        help='Outputfile for the figure and text files')

    #optional parameters
    parser.add_argument('--protname', metavar='protname', type=str,
                        help='Short name for protein of interest.',
                        default="HSA")
    parser.add_argument('--name', metavar='name', type=str,
                        help='Additional name for identification. E.g. Exp1.',
                        default="Exp1")
    parser.add_argument('--annotateK', metavar='annotateK', type=bool,
                        help='Annotate all Lysnes (True) or every 10th residue\
                        (False)',
                        default=False)
    parser.add_argument('--norm', metavar='norm', type=bool,
                        help='Normalize cross-link counts?', default=True)
    parser.add_argument('--residue1', metavar='residue1', type=str,
                        help='cross-linkable residues (1)', default="K")
    parser.add_argument('--residue2', metavar='residue2', type=str,
                        help='cross-linkable residuess (1)', default="STY")


    args = parser.parse_args()
    if (args.in_links is None) or (args.out_path is None):
        print parser.print_usage()
        sys.exit("You need to provide a  link, fasta, and out-file! Check \
your input parameters.")

    else:
        print "Checked input: parameters are fine."
        pass

    #multiple arguments possible
    links_in = args.in_links

    #one input parameter
    fasta_in = args.in_fasta[0]
    name = args.name
    out_path = os.path.abspath(args.out_path[0]) + "/" + name + "_"
    protname = args.protname
    annotateK = args.annotateK
    norm = args.norm
    res1 = args.residue1
    res2 = args.residue2

    print """Input Paramters:
-------------------------------------------"""
    print "Protname: {}".format(protname)
    print "Links (first infile): {}".format(os.path.basename(links_in[0]))
    print "FASTA: {}".format(os.path.basename(fasta_in))
    print ""

    print """Start iteration:
-------------------------------------------"""
    counter = 1
    for links_file in links_in:
        report = CLQC.QC_report("", links_file, fasta_in)
        report.Create_CoverageProfile(out_path, winsize=10, norm=norm,
                                      residue1=res1,
                                      residue2=res2, annotateK=annotateK,
                                      proteinname=protname)
        counter += 1

    print """Report:
-------------------------------------------"""
    print report.report


    print "\r\nFiles written to: \r\n{}".format(out_path+protname+str(counter).zfill(2))




#==============================================================================
# test run
#==============================================================================
"""

ContactMap.py --in_fasta "/home/sven/data/Qutsvs.Velos/1AO6.fasta" --in_links '/home/sven/data/LumosFragmentation2015/CID_5pLink_1.00000_1.00000_1.00000_0.05000_1.00000_100.00000_true_Links_xiFDR1.0.4.csv' --out_path test/ --name HSA
"""

