from os import chdir
import csv

chdir("/home/davis/Desktop/work stuff/final paper/figures")
#chdir(r"C:\Users\Nanda Lab\Desktop\Alex\final paper\figures")

for file_name, moment_name in [('exclusive_moment.csv','ez_exc'),
                               ('inclusive_moment.csv','ez_inc'),
                               ('exc_eisenmoment.csv', 'ei_exc'),
                               ('inc_eisenmoment.csv', 'ei_inc')]:
    try:
        f = open(file_name)
        freader = csv.reader(f)
        for line in freader:
            pdbid = line[0]
            # home version doesn't have 1qd5
            if pdbid == '1QD5':
                continue
            moment = Vector([float(c) for c in line[1:]])
            moments.draw_vector(pdbid + '.' + moment_name,
                                moment, groupdict[pdbid].cored_center)
            groupdict[pdbid].__dict__.update(((moment_name, moment),))
    finally:
        f.close()

print('done')
