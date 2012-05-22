from __future__ import division

class Count(object):
    self.data = None
    def __set__(self, obj, value):
        self.data = value
        none_are_none = True
        aromatic_counts = (obj.tryptophans, obj.phenylalanines,
                           obj.tyrosines, obj.histidines)
        
        for count in aromatic_counts:
            if count is None:
                none_are_none = False
                break

        if none_are_none and obj.residues is not None:
            obj.aromatics = sum(aromatic_counts)
            obj.aromatic_freq = obj.aromatics / obj.residues

     def __get__(self, obj, type=None):
        return self.data


class FrequencyReport(object):
    def __init__(self, pdbid):
        self.pdbid = pdbid

        self.tryptophans = Count()
        self.phenylalanines = Count()
        self.tyrosines = Count()
        self.histidines = Count()
        self.aromatics = Count()

        self.residues = None
        self.aromatic_freq = None

