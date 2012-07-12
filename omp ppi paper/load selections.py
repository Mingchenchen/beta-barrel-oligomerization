# Load selections
# This algorithm probably slows down quickly as you add more objects,
# this is for only a few objects
selection_dir = r'C:\cygwin\home\alex\beta barrels\omp ppi paper\ezbeta coloring\selections'
for pym_obj in cmd.get_names():
    name_match = re.match('aligned_(....)',pym_obj)
    if name_match is None:
        continue
    obj_name = name_match.group(1)

    for selection_file in os.listdir(selection_dir):
        sele_match = re.match('(.*)\.csv', selection_file)
        sele_name = sele_match.group(1)
       
        selection = None
        with open(selection_dir + '/' + selection_file) as f:
            for row in csv.reader(f):
                if row[0].upper() == obj_name.upper():
                    selection = row[1:]
                    break
        if selection == None:
            # It's okay to raise a generic exception here, because
            # this code isn't ina  function so it couldn't possibly have
            # anything catching exceptions
            raise Exception('object name {} not found in selection file {}'\
                            .format(obj_name, selection_file))
        cmd.select(sele_name, 'none')
        for resi in selection:
            cmd.select(sele_name, sele_name + ' | i. ' + resi)
