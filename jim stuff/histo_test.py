from __future__ import division
import csv

class Histogram(list):
    def __init__(self, binsize):
        self.binsize = binsize

    def extend_to(self, i):
        '''x.extend(i): extend so that highest index is i
        All entries added are zeros'''
        if len(self) <=i:
            # Its highest index should be i, so its length should be i+1
            self += [0 for n in range(i+1 - len(self))]

    def __getitem__(self, i):
        '''Retrieve item i, or make list bigger until there is an item i'''
        self.extend_to(i)
        return list.__getitem__(self, i)

        
    def __setitem__(self, i, y):
        '''Set item i, extending list so that item i exists'''
        self.extend_to(i)
        return list.__setitem__(self, i, y)
    
    def add_data(self, value):
        self[value // self.binsize] += 1

    def save(self, filename):
        with open(filename, 'wb') as f:
            csv.writer(f).writerows((i * self.binsize, count) \
                                    for i, count in enumerate(self))
            
