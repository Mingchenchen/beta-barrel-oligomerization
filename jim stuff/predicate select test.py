def select_by_predicate(name, function, arg_list, superset):

    def conditional_add(resi, *args):
        if function(*args):
            cmd.select(name, '{} | ({} & i. {})'\
                             .format(name, superset, resi))
    
    stored.conditional_add = conditional_add

    cmd.select(name, 'none')


    cmd.iterate_state(1, superset,
                      'stored.conditional_add(resi, '\
                      + ', '.join(arg_list) \
                      + ')')



# Select top/bottom halves

start = time.time()
def extracellular_half(z):
    return z>0
def periplasmic_half(z):
    return z<0

for pdbid in groupdict.keys():
    pdbid = pdbid.upper()
    select_by_predicate(pdbid + '.ec_half', extracellular_half,
                       ['z'], pdbid + '.molecule')
    select_by_predicate(pdbid + '.pp_half', periplasmic_half,
                       ['z'], pdbid + '.molecule')
print('top/bottom selections took ' + str(time.time() - start))
