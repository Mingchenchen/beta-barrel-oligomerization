cmd.color('green', '*')

for pdbid in groupdict.keys():
    stored.removed_from_removed = set()
    target = '({0}.exposed & !{0}.cored) & !({0}.removed)'.format(pdbid)
    cmd.color('red', target)
    cmd.iterate(target, 'stored.removed_from_removed.add(resi)')
    print(pdbid + ':')
    print(' '.join(stored.removed_from_removed))
    