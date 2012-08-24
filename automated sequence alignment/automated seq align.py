from Bio.Align.Applications import ClustalwCommandline
import Bio.AlignIO

def call_clustalw(seq_file, mat_file, output_file):
    '''The function that actually produces the sequence alignment.'''
    cline = ClustalwCommandline('clustalw', infile=seq_file,
                                matrix=mat_file, outfile = output_file)

def cluster_retriever(cluster_name):
    '''Returns a multiple sequence alignment corresponding to the given
    cluster name.'''
    if cluster_name[:3].upper() == 'OMP':
        path = 'clusters/{}.clu'.format(cluster_name)
    
    else:
        path = 'clusters/OMP.{}.clu'.format(cluster_name)

    return Bio.AlignIO.read()


