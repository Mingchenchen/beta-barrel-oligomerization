from __future__ import division
from Bio import AlignIO
from sundries import CIDict
from sundries import one_letter
from Bio.PDB import PDBParser
import warnings
import csv

# Retrieve the sequences from the BBTMOUT alignment, including -'s for gaps
bbtm_align= list(AlignIO.read('Swiss-PDB structural alignment.aln',
                          'clustal'))
# Assuming the first is 1A0S, the second is 1AF6:
sequences = CIDict((('1A0S',str(bbtm_align[0].seq)),
                    ('1AF6',str(bbtm_align[1].seq))))

# Check that I'm right about the first being 1AF6, the second being 1A0S
firstfive_of = CIDict()
for pdbid, sequence in sequences.items():
    firstfive_of.update({pdbid: ''})
    for letter in sequence:
        if letter != '-':
            firstfive_of[pdbid] += letter
        if len(firstfive_of[pdbid]) == 5:
            break

assertion_error_message = 'wrong aligned sequences in sequences dictionary'
assertion_error_message += ": 1a0s's first five are {},"\
                           .format(firstfive_of['1a0s'])\
                           +" and 1af6's first five are {}"\
                           .format(firstfive_of['1af6'])

assert firstfive_of['1a0s'] == 'SGFEF' \
       and firstfive_of['1af6'] == 'VDFHG',\
       assertion_error_message

# Retrieve Daniel's aligned structures
structures = CIDict()
parser = PDBParser()
with warnings.catch_warnings():
    # When importing Daniel's aligned structures, the PDBParser gives
    # warnings about "invalid or missing" b factors and occupancies
    # There are so many that if you let it display warnings it'll never
    # finish parsing
    warnings.simplefilter('ignore')
    for pdbid in ('1A0S', '1AF6'):
        structures.update({pdbid: \
                           parser.get_structure(pdbid,
                           'aligned_{}.pdb'.format(pdbid))})
z_coords = CIDict()
for pdbid in structures.keys():
    z_coords.update({pdbid: list()})
    iter_sequence = iter(sequences[pdbid])
    for residue in structures[pdbid].get_residues():
        try:
            calpha = residue.child_dict['CA']
        except KeyError:
            # HOH and other heteroatoms will not have a C-alpha
            # They will also not have a corresponding letter in the sequence
            # So, skip them
            continue

        resi = residue.get_id()[1]
        z = calpha.get_coord()[2]

        # This structure should have the same sequence as the one
        # loaded from the multiple sequence alignment.
        # Raise an error if this residue doesn't match the next letter
        # in the sequence alignment:
        seq_equiv = iter_sequence.next()
        while seq_equiv == '-':
            # Skip gaps:
            seq_equiv = iter_sequence.next()
        assert one_letter[residue.get_resname()] == seq_equiv,\
               'Residue with resi {} in pdb file'.format(resi) \
               + 'doesn\'t match sequence in {}'.format(pdbid)
        
        # Save z coordinate
        z_coords[pdbid].append(z)

z_diff_when_known_is = CIDict()

for known_structure in ('1A0S', '1AF6'):
    if known_structure == '1AF6':
        unknown_structure = '1A0S'
    if known_structure == '1A0S':
        unknown_structure = '1AF6'
    
    print('---')
    print('Known: ' + known_structure)
    print('Unknown: ' + unknown_structure)

    z_diff = list()

    iter_known_z = iter(z_coords[known_structure])
    iter_unknown_z = iter(z_coords[unknown_structure])

    for known_seq_letter, unknown_seq_letter \
        in zip(sequences[known_structure], sequences[unknown_structure]):
        print('')
        print('Pair: '+known_seq_letter+unknown_seq_letter)
        if known_seq_letter == '-':
            # There is no position in the known structure corresponding to
            # this residue in the sequence of the protein of unknown structure.
            # Ignore that residue.
            iter_unknown_z.next()
        elif unknown_seq_letter == '-':
            # A position in the known structure corresponds to a
            # gap. Ignore that position.
            iter_known_z.next()
        else:
            # Predicted position minus real position corresponding to
            # this sequence position
            try:
                predicted = iter_known_z.next() 
                real = iter_unknown_z.next()
                z_diff.append(predicted-real)
                print('predicted {} - real {} = {}'\
                      .format(predicted, real, z_diff[-1]))
                
            except StopIteration:
                continue
        
    z_diff_when_known_is.update({known_structure: z_diff})
    with open('known {}, unknown {} z_diff, swiss pdb structure alignment.csv'\
              .format(known_structure, unknown_structure), 'wb') as f:
        fwriter = csv.writer(f)
        fwriter.writerows([[i] for i in \
                           z_diff_when_known_is[known_structure]])
    

