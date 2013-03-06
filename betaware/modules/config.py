# Copyright (C) 2012 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#  
#  config.py - This file is part of BETAWARE
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

__all__ = ['DESCRIPTION',
           'EPILOG'
           'THRESHOLD',
           'ELM_MODELS',
           'CRF_MODEL_FILE',
           'CRF_DECODING',
           'AA_ORDER',
           'SS_LINE_LEN']


DESCRIPTION = """
 betaware.py: Discrimination and topology prediction of
              trans-membrane beta-barrels with N-to-1 NNs and CRFs.
 
   Copyright (C) 2012 Castrense Savojardo
   Department of Computer Science
   Bologna Biocomputing Group
   University of Bologna, Italy.
   savojard@biocomp.unibo.it
"""

EPILOG = """References: 

[1] Savojardo C., Fariselli P., Casadio R., Improving the detection 
    of transmembrane beta-barrel chains with N-to-1 Extreme Learning Machines, 
    Bioinformatics 27 (22): 3123-3128, 2011.
[2] Fariselli P., Savojardo C., Martelli P.L., Casadio R., 
    Grammatical-Restrained Hidden Conditional Random Fields for 
    Bioinformatics Applications, AlMoB 4:13, 2009.
"""

TH_RNG_MIN = 0.0
TH_RNG_MAX = 0.5

SENS_MIN = 0.0
SENS_MAX = 1.0

ELM_MODELS = [('data/NRPDB.w7.h1600.logistic.model', 7),
              ('data/NRPDB.w11.h1600.logistic.model', 11)]

CRF_MODEL_FILE = 'data/CRFModel.model'

CRF_DECODING = 'posterior-viterbi-sum'

AA_ORDER = 'VLIMFWYGAPSTCHRKQEND'

SS_LINE_LEN = 60

