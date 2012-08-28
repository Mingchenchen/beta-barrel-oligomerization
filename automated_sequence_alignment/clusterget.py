import csv

pdbid_clusterid = list()
with open('structure dataset and clusterguide.csv', 'rb') as f:
    first = True
    for row in csv.reader(f):
        # Skip the first line
        if first:
            first = False
            continue

        # Get the pdbid
        if row[1] != '':
            pdbid = row[1]
        else:
            continue

        # Get the clustername from the Cluster Name column
        if row[2] != '':
            cluster = row[2]
        # But maybe this structure wasn't located to a supercluster
        # In that case, get it from the subcluster column.
        else:
            cluster = row[3]

        pdbid_clusterid.append((pdbid, cluster))

with open('ezbeta clusters/clusters.csv', 'wb') as o:
    csv.writer(o).writerows([[cluster] \
                             for pdbid, cluster in pdbid_clusterid])

