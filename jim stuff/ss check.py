def ss(group, target_resi):
    '''ss(group,resi)
    retrieves secondary structure using non ppi dataset spreadsheet'''
    for resi, ss in group.non_ppi_data.get_columns(['resi','ss']):
        if resi != '' and int(resi)==int(target_resi):
            return ss