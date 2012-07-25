def select_range(bottom, top):
    cmd.select('range', 'none')
    for i in cmd.get_names():
        if ".molecule" not in i:
            continue
        stored.resis = list()
        cmd.iterate(i, "stored.resis.append(resi)")
        resis = stored.resis
        for resi in resis:
            cmd.iterate_state(1, i + ' & n. ca & i. ' + resi,
                              'stored.z = z')
            z = stored.z
            if z > bottom and z < top:
                cmd.select('range',
                           'range | {0} & n. ca & i. {1}'.format(i, resi))
