#chdir(r'C:\Users\Nanda Lab\Desktop\Alex\final paper\figures\magnitudes')
chdir('/home/davis/Desktop/work stuff/final paper/figures/magnitudes')

def savemags(path, moment):
    f = open(path, 'wb')
    fwriter = csv.writer(f)
    
    moment_list = list()
    for name, group in groupdict.items():
        moment_list.append([name, group.__dict__[moment].norm()])
    moment_list = sorted(moment_list, key = lambda x: x[1])[::-1]

    fwriter.writerows(moment_list)
    f.close()

for path, moment in [('ez_exc.csv', 'ez_exc'),
                     ('ez_inc.csv', 'ez_inc'),
                     ('ei_exc.csv', 'ei_exc'),
                     ('ei_inc.csv', 'ei_inc')]:
    savemags(path, moment)

print('done')
