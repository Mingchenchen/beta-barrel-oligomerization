#!/usr/bin/env python
#
# Copyright (C) 2012 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#  
#  betware.py - This file is part of BETAWARE
#  
#  BETAWARE is a free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation; either version 3 of the License,
#  or (at your option) any later version.
# 
#  BETAWARE is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with BioCRF; if not, write to the Free Software Foundation,
#  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.

__version__ = 1.0

import sys
import os

try:
    import argparse
except ImportError:
    sys.stderr.write("Error: Python module argparse not found.\n")
    sys.exit(1)

try:
    import numpy
except ImportError:
    sys.stderr.write("Error: Python module numpy not found.\n")
    sys.exit(1)

try:
    sys.path.append(os.environ['BETAWARE_ROOT'])
except KeyError:
    sys.stderr.write("BETAWARE_ROOT is not set and should point to betaware package root.\n")
    sys.stderr.write("Please, run:\n\n")
    sys.stderr.write("export BETAWARE_ROOT=/betaware/installation/path\n\n")
    sys.stderr.write("and try again. See README for more information.\n")
    sys.exit(1)

import modules.biocompy.slfn as slfn
import modules.biocompy.crf as crf
import modules.utils as utils
import modules.config as config

class BetawareError(Exception):
    """ 
    Generic Betaware exception. 
    """
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def parse_arguments():
    parser = argparse.ArgumentParser(prog = 'betaware.py',
                                     formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description = config.DESCRIPTION,
                                     epilog = config.EPILOG,
                                     usage = '%(prog)s -f FASTA_File -p PROFILE_File [OPTIONS]')
    options = parser.add_argument_group('OPTIONS')


    options.add_argument('-f',
                         help = 'Protein sequence in FASTA format. Required.',
                         type = argparse.FileType('r'),
                         dest = 'fasta',
                         metavar = 'FILE',
                         required = True)
    options.add_argument('-p',
                         help = 'Protein sequence profile. Required.',
                         type = argparse.FileType('r'),
                         dest = 'inputFile',
                         metavar = 'FILE',
                         required = True)
    options.add_argument('-o',
                         help = 'Output prediction file. Optional, default: STDOUT.',
                         type = argparse.FileType('w'),
                         dest = 'outFile',
                         default = sys.stdout,
                         metavar = 'FILE',
                         required = False)
    options.add_argument('-a',
                         help = 'A string specifying the correspondence between amino acids and columns in the sequence profile. Optional, default: %(default)s.',
                         default = config.AA_ORDER,
                         dest = 'aaOrder',
                         metavar = 'STRING')
    options.add_argument('-t',
                         help = 'Always predict topology, also when protein is predicted as non-TMBB. Optional, default = False.',
                         action = 'store_true',
                         dest = 'report_topology')
    options.add_argument('-s',
                         help = 'Sensitivity of the detection algorithm between 0-1. The higher it is the higher is the chance to get false positives. Optional, default: %(default)s.',
                         type = float,
                         default = 0.5,
                         dest = 'sens',
                         metavar = 'VALUE')
    ns = parser.parse_args()
    return ns

def predict_topology(profile):
    crfmodel = crf.CRF()
    crfmodel.parse(os.path.join(os.environ['BETAWARE_ROOT'],
                                config.CRF_MODEL_FILE))
    topo = crfmodel.predict([(profile, None)], algo = config.CRF_DECODING, prob = True)
    labels = [x[0] for x in topo[0]]
    probs = [x[1] for x in topo[0]]
    return ("".join(labels), probs)

def detect_TMBB(profile):
    ensemble = slfn.SLFNEnsemble()
    for model in config.ELM_MODELS:
        ensemble.addModel(os.path.join(os.environ['BETAWARE_ROOT'], model[0]), model[1])
    ensemble.compute_H_row(profile)
    pred = ensemble.run()[0, 0]
    return pred

def main():
    ns = parse_arguments()
    outFile = ns.outFile

    utils.log(config.DESCRIPTION + '\n', outFile)

    profile = utils.read_profile(ns.inputFile, ns.aaOrder, config.AA_ORDER)
    fr = utils.read_single_fasta(ns.fasta)
    if not profile.shape[0] == len(fr):
        raise BetawareError("Error: FASTA sequence length doesn't match profile dimension.")
    pred = detect_TMBB(profile)

    utils.log('Sequence id     : %s\n' % (fr.id,), outFile)
    utils.log('Sequence length : %d\n' % (len(fr.seqdata),), outFile)

    sens = 1.0 - ns.sens
    th = config.TH_RNG_MIN + ((sens - config.SENS_MIN) *
                              (config.TH_RNG_MAX - config.TH_RNG_MIN) /
                              (config.SENS_MAX - config.SENS_MIN))
    if pred >= th:
        utils.log('Predicted TMBB  : Yes\n', outFile)
        topology, probs = predict_topology(profile)
        utils.log('Topology        : %s\n' % (utils.get_TM_segments(topology),), outFile)
        utils.log(utils.get_SS_string(topology, fr.seqdata, probs, l = config.SS_LINE_LEN), outFile)
    else:
        utils.log('Predicted TMBB  : No\n', outFile)
        if ns.report_topology:
            topology, probs = predict_topology(profile)
            utils.log('TMB Strands     : %s\n' % (utils.get_TM_segments(topology, l = 60, o = 18),), outFile)
            utils.log(utils.get_SS_string(topology, fr.seqdata, probs, l = config.SS_LINE_LEN), outFile)
    utils.log('//\n\n', outFile)
    sys.exit(0)

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        raise
        print e
        sys.exit()
