"""All the functions and objects for getting data from csv files and ConSurf
.grades files.

THE GFILE CLASS: initialize with the path of a ConSurf .grades file.
    Can retrieve ConSurf scores, and check whether ConSurf marked
    a residue as having insufficient data for an accurate score
THE SPREADSHEET CLASS:  initialize with the path of a csv file.
    Can retrieve columns.
SPREADSHEET RELATED FUNCTIONS:
    workbook(path, filename_templates, key_template): loads all the csv files
    in a directory as spreadsheets.
    phrasebooks(csvpath): loads a group of phrasebooks from a csv file
MISCELLANEOUS FUNCTIONS:
    file_dict(path, filename_templates, key_template): make a dictionary
        mapping to all files in a given directory that fit given templates"""
    
import csv
import re
import os
    
# THE GFILE CLASS

def gradesplit(item):
    """A function to split column data from a grades file.
    If you just use split(), the confidence interval column gets cut in half whenever its second number is positive.
    Also, I want a separate 'uncertainty flag' column at the end, in addition to the star in the color column."""
    
    output = item.split()
    
    # Merge confidence interval numbers
    if len(output) == 10:
        output[5] += output[6]
        del output[6]
    
    # Add a column at the end with a N indicating an uncertainty flag and an Y indicating none (since the flag indicates uncertainty)
    if '*' in output[4]:
        output.append('N')
    else:
        output.append('Y')
        
    return output

# A dictionary with the indices for values after gradesplit finishes with a line
column = {'res': 2, 'score': 3, 'flag': 9}      
        
class Gfile(object):
    """A consurf .grade file. Initialized with a filename.
    
    get_score_and_flag: Gets tuple with a ConSurf score and uncertainty flag (bool) for a given residue id number
    get_score: Gets a ConSurf score for a given resi
    get_flag:  Gets an uncertainty flag for a given resi
    get_score_if_certain: Returns ConSurf score if certain, 0 if not
    score_list: returns list of tuples, (resi, score). Resis are ints, scores are floats.
    certain_score_list: like scoreList using get_score_if_certain
    flag_list: returns list of tuples, (resi, flag)"""
    # Residue name, number and chain is column['res'], score is column['score'] and confidence flag ('Y' or 'N') is column['flag']
    # NOTETOSELF: If the getters turn out to be too slow, rewrite them with dictionaries
    def __init__(self,filename):
        file = open(filename,'rb')
        try:
            self.data = file.readlines()
            # First 15 lines and last 4 lines are not computable data
            del self.data[:15]
            del self.data[-4:]
            self.data = map(gradesplit,self.data)
            file.close()            
        except:
            self.data = []
            file.close()
    
    def get_score_and_flag(self,resID):
        """Gets a consurf score and flag for a given residue id number, returns as a tuple.
        The flag will return as a boolean, True for confidence, False for insufficient data. -1 indicates an error."""
        output = (-1,-1)
        for i in self.data:
            # Fast check, rules out all but a few possibilities
            if str(resID) in i[column['res']]:          
                # Slow check using regular expressions
                if str(resID) == re.sub("\D", "",i[column['res']]):
                    output = (float(i[column['score']]),i[column['flag']] == 'Y')
                    break
        return output
        
    def get_score(self,resID):
        """Gets a ConSurf score"""
        return self.get_score_and_flag(resID)[0]
        
    def get_flag(self,resID):
        """ Gets an uncertainty flag: True for confidence, False for insufficient data. -1 indicates an error."""
        return self.get_score_and_flag(resID)[1]
        
    def get_score_if_certain(self,resID):
        """Returns score if there's a certainty flag, returns 0 if there's an insufficient data flag"""
        results = self.get_score_and_flag(resID)
        if results[1] is True:
            return results[0]
        else:
            return 0
            
    def score_list(self):
        """Makes a list of tuples, (resID, score)"""
        output = []
        for i in self.data:
            output.append((int(re.sub("\D", "",i[column['res']])), float(i[column['score']])))
        return output
    
    def certain_score_list(self):
        """Makes a list of tuples, (resID, score). If the score is flagged uncertain, it will be recorded as 0."""
        output = []
        for i in self.data:
            # Uncertain score
            if (i[ column['flag'] ]=='N'):
                output.append( ( int(re.sub("\D", "",i[ column['res'] ])), 0 ) )
            # Certain score:
            else:
                output.append( (int(re.sub("\D", "",i[ column['res'] ])), float(i[ column['score'] ])) )
        return output
    
    def flag_list(self):
        """Makes a list of tuples, (resID, uncertainty flag)"""
        output = []
        for i in self.data:
            output.append((int(re.sub("\D", "",i[column['res']])), i[column['flag']] == 'Y'))
        return output
        
    def score_and_flag_list(self):
        """Makes a list of tuples, (resID, score, uncertainty flag)"""
        output = []
        for i in self.data:
            output.append((int('0'+re.sub("\D", "",i[column['res']])), float(i[column['score']]), i[column['flag']] == 'Y'))
        return output

# SPREADSHEET CLASS

def phrasebooks(csvpath):
    """Loads a group of phrasebooks from a file.
    Phrasebook files are csv files. Each phrasebook consists of two columns.
    At the top of the first column is the phrasebook title. Below it is the
    list of keys.
    To the right of each key is the title to which those keys should refer.
    Titles must be in the first line of the csv file.
    ARGUMENTS:
    The path of a csv file.
    RETURNS:
    A dictionary mapping phrasebook titles to phrasebooks."""
    
    output = {}
    csvfile = open(csvpath, 'rb')
    
    try:
        csvreader = csv.reader(csvfile)     
        # Get titles from the first line:
        title_to_column = []
        for index, entry in enumerate(csvreader.next()):
            if entry != '':
                output.update({entry: {}})
                title_to_column.append((index, entry))      
        # Add pairs to the phrasebooks in output:
        for line in csvreader:
            for index, title in title_to_column:
                if line[index] != '':
                    output[title].update({line[index]: line[index+1]})
        csvfile.close()
        
    except:
        csvfile.close()
        raise
        
    return output
        
def translate(mapping, phrasebook):
    """Given a pair of dictionaries, is a generator that can be
    used to create a new dictionary. Using the new
    dictionary is identical to using first_arg[second_arg[key]], where
    first_arg and second_arg are the first and second dictionaries given as
    arguments"""
    for key in phrasebook.iterkeys():
        try:
            yield (key, mapping[phrasebook[key]])
        except KeyError:
            continue
    
class Spreadsheet(object):
    """Initialized with the path of a csv file. Can give a phrasebook as a
    second argument.
    A phrasebook is a dictionary whose items are titles of columns in the
    spreadsheet.
    If a phrasebook is given, the keys of the phrasebook must be used to refer
    to the columns.
    If no phrasebook is given, the column titles must be used to refer to the
    columns.
    get_columns: given a list of column titles, retrieves
        those columns
    get_titles: retrieves mapping between titles and column numbers, good for
        solving mysterious problems
    replace_titles: replaces current titles with those given in a phrasebook.
    add_titles: adds titles in phrasebook to currently known titles.
    replace_title: takes a tuple (newtitle,oldtitle) and replaces oldtitle
        with newtitle
    add_title: takes a tuple (newtitle,oldtitle) and the spreadsheet will
        recognize newtitle as referring to the same column as oldtitle"""
    # Atrributes: self.lines, self.titles
    
    def __init__(self, csvpath, phrasebook = None):
        csvfile = open(csvpath, 'rb')
        try:
            reader = csv.reader(csvfile)
            self.lines = []
            for line in reader:
                self.lines.append(line)
        except:
            print 'fucked up loading ' + csvpath
            csvfile.close()
            raise
        csvfile.close()
        
        # Generate a dictionary that maps column titles to column numbers
        tuples_list = []
        for index, title in enumerate(self.lines[0]):
            tuples_list.append((title, index))
        self.titles = dict(tuples_list)
        
        self.replace_titles(phrasebook)
            
    def get_columns(self, titles_wanted):
        """ Retrieves columns from the spreadsheet.
        ARGUMENT:
        List of column titles.
        RETURNS:
        List of tuples, where tuple n is (nth element of first column title
            given, nth element of second column title given...)"""
        output = []
        for line in self.lines[1:]:
            to_be_appended_list = []
            for title in titles_wanted:
                to_be_appended_list.append(line[self.titles[title]])            
            output.append(tuple(to_be_appended_list))
        return output
        
    def replace_titles(self, phrasebook):
        """Replaces old titles with titles given in a phrasebook.
        A phrasebook is a dictionary that maps the column names to be used
        in the code to the actual column names in the csv file."""
        
        if phrasebook is not None:
            self.titles = dict(translate(self.titles, phrasebook))
            
    def add_titles(self, phrasebook):
        """In addition to current column titles, allows column titles given in
        a phrasebook to be used.
        A phrasebook is a dictionary that maps the column names to be used
        in the code to the actual column names in the csv file."""
        
        if phrasebook is not None:
            self.titles.update((translate(self.titles, phrasebook)))
        
    def get_titles(self):
        """If columns aren't being retrieved as expected, this method
        could help with debugging. Returns the dictionary that maps column
        names to column numbers."""
        return self.titles
        
    def add_title(self, pair):
        """Used to add new column titles without overwriting old ones. So if
        you have 'column_title', you can give this method ('new_title',
        'column_title') and you'll be able to refer to to that column by
        'new_title' or 'column_title'
        ARGUMENT:
        A tuple, with the new title as the first element, and the old title as
        the second.
        RETURNS:
        Nothing"""
        
        # A list containing a single tuple, this has to be the most awkward 
        # possible way to do this:
        self.titles.update([(pair[0], self.titles[pair[1]])])
        
    def replace_title(self,pair):
        """Used to replace a column title with a new one.
        ARGUMENT:
        A tuple, with the new title as the first element, and the old title as
        the second.
        RETURNS:
        Nothing"""
        
        self.add_title(pair)
        del self.titles[pair[1]]

# WORKBOOK SUBROUTINE

def workbook(path, filename_templates, key_template = '(*)'):
    """Creates spreadsheets from all files from a folder and its
    subdirectories whose filenames
    match given templates. Creates a dictionary mapping keys to spreadsheets
    based on a template given for the keys.
    ARGUMENTS:
    path: the folder to be searched
    filename_templates: a list of regular expressions (compiled or as strings)
    key_template: a string. Each instance of (*) is replaced with a group from
        the filename template. The (*)'s are replaced in order, the first one
        being replaced by the first group, etc. I don't like it either, it'd
        be better if you could simply replace them with whatever group you
        like by inserting (*n) where n is a number or something, but I've got
        priorities, I'll add that functionality if I need it. Defaults to
        '(*)'
    RETURNS:
    A dictionary."""
    output = file_dict(path, filename_templates, key_template)
    for key in output.keys():
        output[key] = Spreadsheet(output[key])
    return output

# FILE_DICT SUBROUTINE

def file_dict(path, filename_templates, key_template = '(*)'):
    """Finds file paths for all files in a folder and its subdirectories whose
    names match given templates. Creates a dictionary mapping keys to
    file paths based on a template given for the keys.
    ARGUMENTS:
    path to search
    list of regular expressions to serve as filename templates
    a string with one or more (*)'s in it: this will serve as the template for
        the keys. Each (*) will be replaced with a group from the filename
        templates. Defaults to '(*)'
    RETURNS:
    A dictionary."""
    output = dict()
    for index, string_ in enumerate(filename_templates):
        filename_templates[index] = re.compile(string_)
    for tuple_ in os.walk(path):
        currentdir = tuple_[0]
        for filename in tuple_[2]:
            for rexp in filename_templates:
                match = rexp.match(filename)
                if match is not None:
                    groups = match.groups()
                    key = ''
                    pieces = key_template.split('(*)')
                    for index, segment in enumerate(pieces):
                        key += segment
                        # Possibility for error here if the number of
                        # (*)'s in key_template is more than the number of
                        # .*'s in the filename template.
                        if index != len(pieces) - 1:
                            key += groups[index]
                    output.update({key: currentdir + '/' + filename})
    return output