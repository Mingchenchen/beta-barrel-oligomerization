import csv

# The residue names in the order that I want them:
resns = ['A', 'L', 'I', 'V', 'M', 'D', 'E', 'K', 'N', 'P', 'Q',
         'R', 'H', 'S', 'T', 'F', 'W', 'Y', 'G']

sideways = list()

# Create the first two columns
column1 = ['', 'curve', 'Paq+', 'e0+', 'zmid+', 'n+',
           'Paq-', 'e0-', 'zmid-', 'n-']
sideways.append(column1)
column2 = ['', '', 'Paq+', 'emin+', 'zmin+', 'sigma+',
           'Paq-', 'emin-', 'zmin-', 'sigma-']
sideways.append(column2)


for resn in resns:
    with open('raw param info/ct_{}.csv'.format(resn), 'rb') as f:
        params = list()
        for line in csv.reader(f):
            if line[1] != '':
                params.append(line[1])

        # The first two entries in this column are nTot and nRes, which
        # I don't need
        params = params[2:]
        
        if resn in ['F', 'W', 'Y', 'G']:
            curvetype = 'gaussian'
        else:
            curvetype = 'sigmoidal'

        column = [resn, curvetype] + params
        sideways.append(column)

with open('asymmetric params, bbtmout.csv', 'wb') as f:
    csv.writer(f).writerows(zip(*sideways))
    
        
