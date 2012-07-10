import csv

# The residue names in the order that I want them:
resns = ['A', 'L', 'I', 'V', 'D', 'E', 'K', 'N', 'P', 'Q',
         'R', 'H', 'S', 'F', 'W', 'Y', 'G']

sideways = list()

# Create the first two columns
column1 = ['e0+', 'zmid+', 'n+', 'e0-', 'zmid-', 'n-']
sideways.append(column1)
column2 = ['emin+', 'zmin+', 'sigma-', 'emin-', 'zmin-', 'sigma-']
sideways.append(column2)


for resn in resns:
    with open('raw param info/ct_{}.csv'.format(resn), 'rb') as f:
        params = list()
        for line in csv.reader(f):
            if line[1] != '':
                params.append(line[1])

        # The first two entries in this column are nTot and nRes, which
        # I don't need
        params = params[3:]
        # The first and fifth entries are Paq, which I don't need
        params = params[1:5] + params[6:]
        
        # I'll fill in the curve type (second entry in each column) manually
        column = [resn, ''] + params
        sideways.append(column)

with open('asymmetric params, bbtmout.csv', 'wb'):
    csv.writer(f).writerows(zip(*sideways))
    
        
