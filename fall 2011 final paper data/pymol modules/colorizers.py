from biodata import *
import functools
import color_b
from pymol import stored
from pymol import cmd

def consurf_color(sele, grades, group, gradient = 'red_white_blue',
                  insufficient_data_color = 'yellow', nbins = 30):

    
        group.bsource = grades
      
        # score_and_flag_list returns tuples whose third element is a boolean
        # (the 'flag'). False means that consurf considers there to be
        # insufficient data to decide upon a conservation score.
        list_ = filter(lambda x: x[2], grades.score_and_flag_list())
        list_ = [(x[0],x[1]) for x in list_]

        
        colorize(sele, list_, group,
                 gradient = gradient, no_data_color = insufficient_data_color,
                 nbins = nbins, bsource = grades)
            
def add_consurf_color(group, sele, grades):
    """Adds a consurf_scores attribute and a consurf_color() method to a group
    object. This method will colorize groupname.molecule and groupname.pdb.
    consurf_color has keyword arguments "spectrum" and
    "insufficient_data_color". Spectrums are anything Pymol's spectrum command
    can use, the insufficient_data_color is any color Pymol recognizes.
    Note that unlike the colorizer functions that work from spreadsheets, this
    adds the gfile to the group. This is because multiple colorizations can be
    made from a single spreadsheet, but I can't imagine any need to use one
    Gfile for different colorizations of the same object.
    
    ARGUMENTS:
    Group object
    Gfile object
    
    RETURNS:
    Nothing"""
    
    #NOTETOSELF:
    """    There's an optional moleculenames keyword argument to consurf_color; by
    default it's a list ['moleculenames', 'pdb']. So if you want to color
    'groupname.extra', just call group.consurf_color(moleculenames=['extra'])"""
    # Put that feature back in at some point, and add it to ezbeta moments too
    group.consurf_scores = grades
    
    group.consurf_color = functools.partial(consurf_color,
                                            sele, group.consurf_scores,
                                            group = group)
    
def ezbeta_color(sele, spreadsheet, group, gradient='red_white_blue',
                 no_data_color='yellow', nbins=30):
    colorize(sele, spreadsheet.get_columns(['resi','ez']), group,
             gradient=gradient, no_data_color=no_data_color, nbins=nbins,
             bsource = spreadsheet)
             
def add_ezbeta_color(group, sele, spreadsheet):
    """ See documentation for add_consurf_color."""
    group.ezbeta_color = functools.partial(ezbeta_color, sele, spreadsheet,
                                           group)
            
def colorize(sele, pairslist, group, gradient = 'red_white_blue',
             no_data_color = 'yellow',
             nbins = 30, in_python = False, bsource = None):
    
    """Colorizes given selection based on a list of tuples (resi, number).
    
    ARGUMENTS:
    sele, pairslist, gradient = 'red_white_blue', no_data_color = 'yellow', nbins = 30
    
    RETURNS:
    Nothing"""
    
    # This makes me feel better whenever I'm gonna be appending things to
    # selection strings:
    sele = '(' + sele + ')'

    # Save some time on updating b values by not bothering if it already has
    # the right b values
    if bsource is None or 'bsource' not in dir(group)\
       or group.bsource is not bsource:
        changeb = True
    else:
        changeb = False
    
    stored.in_pairslist = {}
    cmd.iterate(sele, 'stored.in_pairslist.update({resi: False})')
    for resi, number in pairslist:
        if changeb:
            cmd.alter(sele + '& i. ' + str(resi), 'b = ' + str(number))
        # KeyError here if a residue in the pairslist isn't in the molecule
        stored.in_pairslist[str(resi)] = True
    for resi, is_in_pairslist in stored.in_pairslist.iteritems():
        if not is_in_pairslist:
            cmd.color(no_data_color, sele + ' & i. ' + resi)
            sele += ' & (! i. ' + resi + ')'

    if bsource is not None:
        group.bsource = bsource
                
    cmd.spectrum('b',gradient,sele)       
    
def color_selected(part, whole, partcolor, wholecolor):
    pass

#groupdict = groups_from_folder(r'C:\Users\Nanda Lab\Desktop\Alex\consurf',['aligned_(.*).pdb'])
#load_spreadsheets(groupdict, r'C:\Users\Nanda Lab\Desktop\Alex\Excel data', ['(.*)_DeltaBurial.csv'], phrasebook = phrasebooks(r'C:\Users\Nanda Lab\Desktop\Alex\scripts\phrasebooks.csv')['dan_ezbeta'], attribute_name = 'ezbeta_spreadsheet')  
#add_ezbeta_color(groupdict['1A0S'], '1a0s.molecule', groupdict['1A0S'].ezbeta_spreadsheet)
#groupdict['1A0S'].ezbeta_color()
