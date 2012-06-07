import os
import re
import csv

def reduce(source, target):
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
        reader = csv.reader(t)
        writer = csv.writer(t)

        # Find which column contains the residue identifiers:
        title_line = reader.next()

