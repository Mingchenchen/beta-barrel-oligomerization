
test = dict()
for i in ['A', 'B']:
    test.update({i: dict((j, None) for j in ['A', 'B'])})


to_mat(test, ['A', 'B'])

# break

test['A']['A'] = 10
test['A']['B'] = 100
test['B']['A'] = 1000
test['B']['B'] = 10000
to_mat(test, ['A', 'B'])
pi = {'A': 10, 'B': 1}
lam = 1/math.log10(math.e)
b = scoring_matrix(test, lam, pi)

# break for output
