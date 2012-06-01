import biodata
import sundries

# Get the phrasebook for interacting with Daniel's "weights" files
phrasebooks = biodata.phrasebooks('weights phrasebook.csv')

# Retrieve names of csv files containing the resis
# is there something I can use that's simpler than os.walk,
# and just returns all the filenames in a given directory?
weight_file_filenames

# Make spreadsheets
spreadsheets = [biodata.Spreadsheet(filename,
                                    phrasebook = phrasebooks['weights'])\
                for filename in weight_file_filenames]

# weights maps pdbids to spreadsheets
weights = CIDict()
for spreadsheet in spreadsheets:
    # I expect there to only be one pdbid
    pdbids = spreadsheet.get_column

# selections maps pdbids to sets of resis
selections = sundries.CIDict()
