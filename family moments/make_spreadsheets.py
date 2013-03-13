print(x)
for pdbid, family in families(params='normal random parameters.csv').items():
    sel = Selection(family, lambda x: x.rel_acc > .2 and x.ss == 'E')
    dir_ = 'random parameters, relacc over .2, beta sheet, gonnet aligned'
    sel.spreadsheet('{}/{}.csv'.format(dir_, pdbid))
