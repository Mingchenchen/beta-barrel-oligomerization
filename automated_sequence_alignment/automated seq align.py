from alignments import *


def align_all(matrix, output_dir):
    # Retrieve pdbids of structures in the dataset
    # and their associated clusters, from the information in Daniel's
    # thesis.
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

    # Make the alignments
    for pdbid, cluster in pdbid_clusterid:
        # align(output_name, cluster, matrix, *pdbpaths)
        align('{}/{} with {}.clu'.format(output_dir, pdbid, cluster),
              cluster, 'gonnet',
              'ezbeta aligned structures/aligned_{}.pdb'.format(pdbid))
