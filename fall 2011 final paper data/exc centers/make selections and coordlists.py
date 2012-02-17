from os import chdir
import csv

# run /home/davis/Desktop/work stuff/final paper/figures/exc centers/make selections and coordlists.py
chdir('/home/davis/Desktop/work stuff/final paper/figures/exc centers')

def get_exp_sheets():
    stored.exp_sheets = dict()
    for name in groupdict.keys():
        stored.exp_sheets.update(((name, list()),))
        cmd.iterate(name + '.exp_sheets & n. CA',
                    'stored.exp_sheets["%s"].append(resi)' % name)
    return stored.exp_sheets

def take_intersection(selection):
    try:
        f = open('cored 1.csv', 'rb')
        reader = csv.reader(f)
        intersection = dict()
        for line in reader:
            name = line[0]
            intersection.update(((name, list()),))
            for resi in line[1:]:
                if resi in selection[name]:
                    intersection[name].append(resi)
    finally:
        f.close()

    return intersection
    
def make_coord_lists(selections, folder = 'coordlists'):
    for name, resis in selections.items():
        stored.coordlist = list()
        for resi in resis:
            cmd.iterate_state(1, '%s.molecule' % name +\
                              '& n. ca & i. ' + resi,
                              'stored.coordlist.append((x,y))')
        try:
            f = open('%s/%s.csv' % (folder, name), 'wb')
            writer = csv.writer(f)
            writer.writerows(stored.coordlist)
        finally:
            f.close()        
        print('done with ' + name)

    print('done')

def save_selections(selections, file_name = 'exc sheets selections.csv'):
    try:
        f = open(file_name, 'wb')
        writer = csv.writer(f)
        writer.writerows([[name] + resis for name, resis in \
                                             selections.items()])
    finally:
        f.close()

def load_center(file_name, center_name):
    try:
        f = open(file_name, 'rb')
        for line in csv.reader(f):
            groupdict[line[0]].__dict__.update({center_name:
                                                Vector([float(x) \
                                                for x in line[1:]])})
    finally:
        f.close()

def test_centers(selection, center_name, single = True):
    for name, group in groupdict.items():
        stored.sum_ = Vector((0,0,0))
        for resi in selection[name]:
            cmd.iterate_state(1,
                          'n. CA & %s.molecule & i. %s' % (name, str(resi)),
                          'stored.pos = Vector((x,y,z)) ' +
                                      ' - groupdict["%s"].' % name +\
                                      center_name)
            stored.pos = Vector([stored.pos[0], stored.pos[1], 0])
            stored.sum_ += stored.pos.normalize()
        print('%s: %s' % (name, str(stored.sum_)))
        if single:
            break

#load_center('exc centers.csv', 'exc_center')
test_centers(intersection, 'exc_center', single=False)
