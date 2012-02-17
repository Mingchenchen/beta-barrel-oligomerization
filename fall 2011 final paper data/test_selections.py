def test_selections(path):
    try:
        f = open(path)
        reader = csv.reader(f)
        for line in reader:
            name = line[0]
            for resi in line[1:]:
                cmd.color('red', name + '.molecule & i. ' + resi)
    finally:
        f.close()
            
test_selections(r'C:\Users\Nanda Lab\Desktop\Alex\final paper\beta_selections.csv')