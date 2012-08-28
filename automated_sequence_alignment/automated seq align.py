from Bio.Align.Applications import ClustalwCommandline
import Bio.AlignIO
import Bio.SeqIO
import Bio.PDB
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import warnings
import os
import copy
import csv

from sundries import one_letter

def call_clustalw(seq_file, mat_file, output_file):
    '''The function that actually produces the sequence alignment.'''
    cline = ClustalwCommandline('clustalw', infile=seq_file,
                                matrix=mat_file, outfile = output_file)
    
    return cline()

def cluster_retriever(cluster_name):
    '''Returns a multiple sequence alignment corresponding to the given
    cluster name.'''
    if cluster_name[:3].upper() == 'OMP' or 'cluster' in cluster_name:
        path = 'clusters/{}.clu'.format(cluster_name)
    
    else:
        path = 'clusters/OMP.{}.clu'.format(cluster_name)

    return Bio.AlignIO.read(path, 'clustal')

def sequence_retriever(pdb_paths):
    # Match a pathname, with a group that includes the last slash or 
    # backslash and everything after it
    file_from_path = re.compile(r'^.*([/\\].*?)$')

    output = list()

    for path in pdb_paths:
        # Come up with a name to identify the sequence by
        match = re.match(file_from_path, path)
        if match is None:
            name = path
        else:
            name = match.group(0)
        # I think ClustalX ignores everything after a space which is
        # annoying, so replace the spaces with underscores
        name = name.replace(" ", "_")
            
        # Load the structure. See what happens if you try to load one
        # of Daniel's aligned structures without hiding warnings. It's not
        # even about anything significant, for this kind of thing, 
        # just missing b factors and occupancies
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            stru = Bio.PDB.PDBParser().get_structure(name, path)
            
        # Get the sequence from the structure
        seq_as_str = str('')
        for r in stru.get_residues():
            # one_letter is a case-insensitive dictionary mapping
            # 3-letter residue names to 1-letter residue names
            # Biopython PDB structures use 3-letter residue names
            standard_20 = [n.upper() for n in one_letter.keys()]
            resn = r.get_resname().upper()

            # Only take stuff that's actually a residue, and not like
            # some bound lipid or something
            if resn in standard_20:
                seq_as_str += one_letter[resn]

            # In the PDB entry for 1FEP, the FASTA file just has "M"
            # where there's selenomethionine in the structure, 
            # which makes sense because Vik
            # says that putting in the selenium is something that's done
            # to make structure determination easier, a bacterium's FepA
            # is not really going to be full of selenium.
            # So, treat selenomethionine as methionine:
            elif resn == 'MSE':
                seq_as_str += 'M'
        
        seq_as_Seq = Seq(seq_as_str, alphabet = IUPACProtein())

        # Append a SeqRecord to the output
        output.append(SeqRecord(seq_as_Seq, name=name, id=name,
                      description='sequence taken from a structure file'))
        
    return output

def combine(output_filename, msa, seqs):
    output_seqs = list(seqs)

    for record in msa:
        # An HHOMP alignment has lots of useful stuff besides sequences,
        # which I want to ignore here
        # I don't know what "gi" means but all the real sequences seem to
        # have it at the beginning of their id, so,
        # Filter for only those sequences whose id's begin with "gi"
        if record.id[:2] == 'gi':
            # Get rid of gaps
            gapless_as_str = ''.join(filter(lambda x: x!= '-', record))
            gapless_seq = Seq(gapless_as_str)
            gapless_record = copy.deepcopy(record)
            gapless_record.seq = gapless_seq
            output_seqs.append(gapless_record)
    
    Bio.SeqIO.write(output_seqs, output_filename, 'fasta')

def bbtm_filename(t):
    # Make string of t padded with 0's:
    t_string = (3-len(str(t)))*'0' + str(t)

    # Create filename of one of the matrices David Jimenez-Morales
    # sent me
    filename = os.getcwd() + r'\TMout\bbtmout\bbTMout'+t_string

    return filename


def write_series(id_t_triplets, output_filename):
    # Information on the series file format can be found here:
    # http://bips.u-strasbg.fr/en/Documentation/ClustalX/
    # In case the link breaks, the Wayback Machine has an archived page
    # from July 27 2011
    with open(output_filename, 'w') as o:
        o.write('CLUSTAL_SERIES\n')
        o.write('\n')
        for id_bottom, id_top, t in id_t_triplets:
            o.write(' '.join(["MATRIX", str(id_bottom), str(id_top),
                              bbtm_filename(t)]))
            o.write('\n')

# The API
class UnrecognizedMatrix(Exception):
    pass

def align(output_name, cluster, matrix, *pdbpaths):
    # Match a pathname, with a group that includes everything before the 
    # last slash or backslash
    dir_from_path = re.compile(r'^(.*)[/\\].*?$')
    dir_match = re.match(dir_from_path, output_name) 
    if dir_match is not None:
        target_dir = dir_match.group(1)
    else:
        target_dir = ''

    # Set infile parameter
    infile_param = target_dir + '/seqs.fasta'
    cluster_msa = cluster_retriever(cluster)
    pdb_sequences = sequence_retriever(pdbpaths)
    combine(infile_param, cluster_msa, pdb_sequences)
    
    # Set matrix parameter
    # It could be a matrix built into clustal:
    if matrix.lower() in ('blosum', 'pam', 'id', 'gonnet'):
        matrix_param = matrix
    # Or it could be BBTM:
    elif matrix.lower() == 'bbtm':
        # These pairings of sequence identities and evolutionary distances
        # are taken from a 1994 Clustal paper. They're the same as the
        # pairings in that paper 
        # of sequence identities and PAM numbers for the PAM matrices
        series_path = target_dir + r'\bbtm_series.txt'
        write_series([( 0, 40 , 350),
                      (41, 60 , 120),
                      (61, 80 , 60 ),
                      (81, 100, 20 )], series_path)
        matrix_param = os.getcwd() + '\\' + series_path
    # But otherwise I don't know what it is
    else:
        raise UnrecognizedMatrix('Matrix "{0}" not '.format(matrix) \
                                 + 'recognized. Try BBTM,  PAM, '\
                                 + 'GONNET, BLOSUM or ID.')

    return call_clustalw(infile_param, matrix_param, output_name)


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
