import matrices

# Published bbtm ordering including doubtful amino acids and 
# * column:
complete_bbtm_ordering = matrices.published_bbtm_ordering \
                         + ['B', 'Z', 'X', '*']

def write_matrix(mat, ordering, output_filename, title,
                 comments = None, round_ = False):
    with open(output_filename, 'w') as o:
        # Write the first comment line, containing the title
        o.write('# ' + title)

        # Write the comments
        if comments != None:
            for line in comments.split('\n'):
                o.write('# ' + line + '\n')
        
        # Write the row names
        o.write(' '*4 + (' '*3).join(ordering) + '\n')

        # Write the actual entries of the matrix
        for row in ordering:
            o.write(row)
            for col in ordering:
                # Round the entry, if round_ is True
                if round_:
                    entry_as_num = matrices.\
                                   round_to_nearest_int(mat[row,col])
                else:
                    entry_as_num = mat[row,col]
                
                # Create a string that is the entry padded with spaces
                # on the left so that it is four characters long
                entry = str(entry_as_num)
                padded_entry = ' '*(4-len(entry)) + entry
                o.write(padded_entry)
            o.write('\n')

def remake_with_star(in_path, out_path, title, comments='* added'):
    '''Parses one of David Jimenez-Morales's BBTM matrices. Adds a 
    * column, which is equal to the minimum entry in the matrix.
    '''
    # Retrieve the matrix as a dictionary
    mat = matrices.parse(in_path)

    # Find the minimum entry
    min_entry = min(mat.values())

    # Make keys (i, *) and (*, i) map to that minimum entry for all
    # letters i in complete_bbtm_ordering (including (*, *))
    for i in complete_bbtm_ordering:
        mat.update({(i, '*'): min_entry})
        mat.update({('*', i): min_entry})

    # Write the matrix
    write_matrix(mat, out_path, title=title, comments=comments)
