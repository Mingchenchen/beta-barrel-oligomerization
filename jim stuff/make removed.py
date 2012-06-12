for pdbid, group in groupdict.items():
    pdbid = pdbid.upper()
    cmd.select(pdbid + '.removed', '{0}.exposed & ! {0}.cored'.format(pdbid))