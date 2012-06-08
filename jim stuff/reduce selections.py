import os
import re
import csv

def reduce_sheet(source, target):
    '''Daniel's spreadsheets contain the residue names of the protein,
    as well as all the homology models. This makes them large and slow to
    work with.
    I never interact with the residue names from these files, so all the
    information I actually use is simply repeated over and over again.
    This function takes one of Daniel's spreadsheets, and saves a reduced
    form that contains each residue only once, by copying the spreadsheet
    until it finds the first line whose residue id is LOWER than that of
    the previous line; at that point, it ceases the transfer.
    Function of two filenames: reduce(source, target)'''
    
    with open(source, 'rb') as s, open(target, 'wb') as t:
        reader = csv.reader(s)
        writer = csv.writer(t)

        # Find which column contains the residue identifiers:
        title_line = reader.next()
        resi_col = title_line.index('ResIndex')
        last_resi = 0
        for line in reader:
            resi = int(line[resi_col])
            if resi < last_resi:
                print('last resi in {} is {}'.format(line[0], last_resi))
                break
            writer.writerow(line)
            last_resi = resi

for dataset in ('ppi', 'non ppi'):
    for filename in os.listdir(dataset + ' residues'):
        match = re.match('(\d...)\.csv', filename)
        if match is not None:
            reduce_sheet('{} residues/{}.csv'\
                         .format(dataset, match.group(1)),
                         '{} residues reduced/{}.csv'\
                         .format(dataset, match.group(1)))

