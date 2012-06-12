for pdbid in groupdict.keys():
    pdbid = pdbid.upper()
    cmd.select(pdbid + '.not_removed',
               '{0}.molecule & !{0}.removed'.format(pdbid))