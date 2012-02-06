import csv

from sundries import one_letter
stored.one_letter = one_letter


output = list()

for name, group in groupdict.items():
    stored.sequence = ''
    cmd.iterate('{}.molecule & n. ca & polymer'.format(name),
                'stored.sequence += stored.one_letter[resn]')
    output.append((name, stored.sequence))


with open(r'C:\cygwin\home\alex\temp\sequences.csv', 'wb') as f:
    fwriter = csv.writer(f)
    fwriter.writerows(output)
