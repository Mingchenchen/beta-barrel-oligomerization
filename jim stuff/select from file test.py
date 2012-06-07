import csv
import time

start = time.time()

# A faster non_ppi selector that loads it from a file

# Select residues in the non_ppi dataset

# Open the selection from  csv file that is in the format
# pdbid, resi1, resi2, resi3... resiN
non_ppi_selection = CIDict()
with open('selections/non ppi.csv', 'rb') as f:
    for line in csv.reader(f):
        non_ppi_selection.update({line[0]: line[1:]})
# Carry out the selection
for group in groupdict.values():
    sele_name = group.name + '.non_ppi'
    cmd.select(sele_name, 'none')
    for resi in non_ppi_selection[group.name]:
        sele_string = '{} | {} & i. {}'.format(sele_name, group.name, resi)
        try:
            cmd.select(sele_name, sele_string)
        except CmdException:
            print('fucked up trying selection string ' + sele_string)
            raise

print('took ' + str(time.time() - start))
