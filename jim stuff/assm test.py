import csv


with open('non ppi residues/1A0S.csv', 'rb') as s:
    reader = csv.reader(s)

    # Find which column contains the residue identifiers:
    title_line = reader.next()
    resi_col = title_line.index('ResIndex')

    # Make a list of the size of each "segment", where a segment is a
    # sequence of lines where resi increases from the last in each line
    seg_sizes = list()

    last_resi = 0
    cur_seg_size = 0

    for line in reader:
        resi = int(line[resi_col])
        cur_seg_size += 1
        if resi < last_resi:
            seg_sizes.append(cur_seg_size)
            cur_seg_size = 0
            last_resi = 0
        last_resi = resi
    
    seg_sizes.append(cur_seg_size)

for i in seg_sizes:
    print(i)
