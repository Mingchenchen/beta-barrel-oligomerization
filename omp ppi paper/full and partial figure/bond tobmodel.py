stored.resis = list()
cmd.iterate('*', 'stored.resis.append(int(resi))')
for i in stored.resis:
    if i+1 in stored.resis:
        cmd.bond('i. ' + str(i), 'i. '+ str(i+1))
