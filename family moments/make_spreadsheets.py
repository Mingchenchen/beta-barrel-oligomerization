for pdbid, family in families().items():
    sel = Selection(family, lambda x: x.rel_acc > .2 and x.ss == 'E')
    dir_ = 'relacc over .2, beta sheet, gonnet aligned'
    sel.spreadsheet('{}/{}.csv'.format(dir_, pdbid))
