from sundries import CIDict

def selections_by_index(iterable):
    selections = CIDict()
    test_sequences = CIDict()
    for line in iterable:
        selections.update({line[0]: [int(x) for x in line[2:]]})
        test_sequences.update({line[0]: line[1]})
    return selections, test_sequences

def selections_by_resi(iterable):
    selections = CIDict()
    for line in iterable:
        selections.update({line[0]: [int(x) for x in line[1:]]})
    return selections