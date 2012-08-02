with open('parser test bare.txt', 'r') as i,\
     open('parser test for mathematica.txt', 'w') as o:
    o.write('{\n')
    output_lines = list()
    for line in i:
        output_lines.append('{' + ','.join(line.split()) + '}')
    o.write(',\n'.join(output_lines))
    o.write('\n}\n')
