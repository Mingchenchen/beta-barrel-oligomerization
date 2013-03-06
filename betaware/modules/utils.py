# Copyright (C) 2012 Castrense Savojardo
#     Bologna Biocomputing Group
#     University of Bologna, Italy
#     savojard@biocomp.unibo.it
#  
#  utils.py - This file is part of BETAWARE
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

__all__ = ['FastaRec',
           'compute_profile_from_SS',
           'read_profile',
           'read_single_fasta',
           'get_TM_segments',
           'get_SS_string',
           'log']

import numpy

class FastaRec(object):
    def __init__(self, id, seqdata = ''):
        self.id = id
        self.seqdata = seqdata

    def __len__(self):
        return len(self.seqdata)

    def append_seq(self, seqdata):
        self.seqdata.append(seqdata)


def compute_profile_from_SS(sequence, aaOrder):
    l = len(sequence)
    p = numpy.zeros((l, len(aaOrder)))
    for i in range(len(sequence)):
        j = aaOrder.index(sequence[i])
        p[i, j] = 1.0
    return p

def read_profile(handle, aaOrder, outAAOrder):
    lines = handle.readlines()
    l = len(lines)
    p = []
    pos = [aaOrder.index(a) for a in outAAOrder]
    for line in lines:
        line = line.split()
        if len(line) > 0:
            p.append([float(line[i].rstrip('\n')) for i in pos])
    return numpy.array(p)

def read_single_fasta(handle):
    seqres = False
    seq = ''
    id = None
    for line in handle.readlines():
        if line[0] == '>':
            if seqres:
                break
            else:
                id = line[1:].rstrip('\n')
                seqres = True
        else:
           seq += line.rstrip('\n')
    fr = FastaRec(id, seq)
    return fr

def get_TM_segments(sstring, l = 80, o = 20):
    import re
    segs = ''
    al = 0
    for i in re.finditer('T+', sstring):
        if al + len(str((i.start() + 1, i.end()))) > l:
            segs += '\n'
            segs += ' ' * o
            al = 0
        segs += "%d-%d," % (i.start() + 1, i.end())
        al += len(str((i.start() + 1, i.end())))
    return segs[:-1]

def get_SS_string(sstring, seq, probs, l = 80):
    i = 0
    ret = ''
    rngprob = [0.1 * (j + 1) for j in range(0, 10)]
    rngalph = ['j', 'i', 'h', 'g', 'f', 'e', 'd', 'c', 'b', 'a']

    while i < len(seq):
        j = min(i + l, len(seq))
        ret += 'Seq : ' + seq[i:j] + '\n'
        ret += 'SS  : ' + sstring[i:j] + '\n'
        found = False
        probstr = ''
        for pr in probs[i:j]:
            if pr > 1.0:
                pr = 1.0
            if pr < 0.0:
                pr = 0.0
            for k in range(0, 10):
                if pr <= rngprob[k]:
                    probstr += rngalph[k]
                    break
        ret += 'Prob: ' + probstr + '\n'
        ret += '-' * (6 + min(l, len(seq) - i)) + '\n'
        i = j
    return ret

def log(str, handle, required_verbose = 0, verbose = 0):
    if verbose >= required_verbose:
        handle.write(str)
        handle.flush()




