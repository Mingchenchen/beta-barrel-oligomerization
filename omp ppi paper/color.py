import csv

def set_b(filename):
    cmd.select('on_blist', 'none')
    with open(filename, 'rb') as f:
        for resi, b in csv.reader(f):
            b = abs(float(b))
            cmd.select('on_blist', 'on_blist | resi ' + resi)           
            cmd.alter('resi ' + resi,
                      'b = ' + str(b))

def color_b(n=2):
    cmd.spectrum('b','blue_white_red', minimum=0, maximum=n)
    cmd.color('black', '!on_blist')
