def select_by_functions(names, functions, arg_list, superset,
                        group):
    for name, function in zip(names, functions):

        def conditional_add(resi, *args):
            if function(*args):
                cmd.select(name, '{} | ({} & i. {})'\
                                 .format(name, superset, resi))
        
        stored.conditional_add = conditional_add

        cmd.select(name, 'none')


        cmd.iterate_state(1, superset,
                          'conditional_add(resi, '\
                          + ', '.join(arg_list) \
                          + ')')


def ss(group, target_resi):
    '''ss(group,resi)
    retrieves secondary structure using non ppi dataset spreadsheet'''
    for resi, ss in group.non_ppi.get_columns(['resi','ss']):
        if resi != '' and int(resi)==int(target_resi):
            return ss


