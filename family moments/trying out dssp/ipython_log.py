# IPython log file

get_ipython().system(u'dir /on ')
import subprocess
x = subprocess.output(['dssp', '1A0S.pdb'])
x = subprocess.check_output(['dssp', '1A0S.pdb'])
x = subprocess.check_output(['dssp', '1A0S.pdb', '1a0s home calculated.dssp'])
x
#[Out]# ''
import urllib
help(urllib)
dir(urllib)
#[Out]# ['ContentTooShortError',
#[Out]#  'FancyURLopener',
#[Out]#  'MAXFTPCACHE',
#[Out]#  'URLopener',
#[Out]#  '__all__',
#[Out]#  '__builtins__',
#[Out]#  '__doc__',
#[Out]#  '__file__',
#[Out]#  '__name__',
#[Out]#  '__package__',
#[Out]#  '__version__',
#[Out]#  '_ftperrors',
#[Out]#  '_have_ssl',
#[Out]#  '_hexdig',
#[Out]#  '_hextochr',
#[Out]#  '_hostprog',
#[Out]#  '_is_unicode',
#[Out]#  '_localhost',
#[Out]#  '_noheaders',
#[Out]#  '_nportprog',
#[Out]#  '_passwdprog',
#[Out]#  '_portprog',
#[Out]#  '_queryprog',
#[Out]#  '_safe_map',
#[Out]#  '_safe_quoters',
#[Out]#  '_tagprog',
#[Out]#  '_thishost',
#[Out]#  '_typeprog',
#[Out]#  '_urlopener',
#[Out]#  '_userprog',
#[Out]#  '_valueprog',
#[Out]#  'addbase',
#[Out]#  'addclosehook',
#[Out]#  'addinfo',
#[Out]#  'addinfourl',
#[Out]#  'always_safe',
#[Out]#  'basejoin',
#[Out]#  'c',
#[Out]#  'ftpcache',
#[Out]#  'ftperrors',
#[Out]#  'ftpwrapper',
#[Out]#  'getproxies',
#[Out]#  'getproxies_environment',
#[Out]#  'getproxies_registry',
#[Out]#  'i',
#[Out]#  'localhost',
#[Out]#  'main',
#[Out]#  'noheaders',
#[Out]#  'os',
#[Out]#  'pathname2url',
#[Out]#  'proxy_bypass',
#[Out]#  'proxy_bypass_environment',
#[Out]#  'proxy_bypass_registry',
#[Out]#  'quote',
#[Out]#  'quote_plus',
#[Out]#  'reporthook',
#[Out]#  'socket',
#[Out]#  'splitattr',
#[Out]#  'splithost',
#[Out]#  'splitnport',
#[Out]#  'splitpasswd',
#[Out]#  'splitport',
#[Out]#  'splitquery',
#[Out]#  'splittag',
#[Out]#  'splittype',
#[Out]#  'splituser',
#[Out]#  'splitvalue',
#[Out]#  'ssl',
#[Out]#  'string',
#[Out]#  'sys',
#[Out]#  'test',
#[Out]#  'test1',
#[Out]#  'thishost',
#[Out]#  'time',
#[Out]#  'toBytes',
#[Out]#  'unquote',
#[Out]#  'unquote_plus',
#[Out]#  'unwrap',
#[Out]#  'url2pathname',
#[Out]#  'urlcleanup',
#[Out]#  'urlencode',
#[Out]#  'urlopen',
#[Out]#  'urlretrieve']
get_ipython().magic(u'pinfo urllib.urlretrieve')
urlretrieve('ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/1a0s.dssp', filename="1a0s from dssp database.dssp")
urllib.urlretrieve('ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/1a0s.dssp', filename="1a0s from dssp database.dssp")
#[Out]# ('1a0s from dssp database.dssp', <mimetools.Message instance at 0x036DD878>)
home  = open('1a0s home calculated.dssp').readlines()
home[:10]
#[Out]# ['==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000 ==== DATE=15-FEB-2013     .\n',
#[Out]#  'REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .\n',
#[Out]#  'HEADER    OUTER MEMBRANE PROTEIN                  07-DEC-97   1A0S                                                             .\n',
#[Out]#  'COMPND   2 MOLECULE: SUCROSE-SPECIFIC PORIN;                                                                                   .\n',
#[Out]#  'SOURCE   2 ORGANISM_SCIENTIFIC: SALMONELLA TYPHIMURIUM;                                                                        .\n',
#[Out]#  'AUTHOR    K.DIEDERICHS,W.WELTE                                                                                                 .\n',
#[Out]#  ' 1239  3  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .\n',
#[Out]#  ' 44143.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n',
#[Out]#  '  910 73.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .\n',
#[Out]#  '    6  0.5   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n']
database = open('1a0s from dssp database.dssp').readlines()
database[:10]
#[Out]# ['==== Secondary Structure Definition by the program DSSP, CMBI version by M.L. Hekkelman/2010-10-21 ==== DATE=2012-03-20        .\n',
#[Out]#  'REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .\n',
#[Out]#  'HEADER    OUTER MEMBRANE PROTEIN                  07-DEC-97   1A0S                                                             .\n',
#[Out]#  'COMPND   2 MOLECULE: SUCROSE-SPECIFIC PORIN;                                                                                   .\n',
#[Out]#  'SOURCE   2 ORGANISM_SCIENTIFIC: SALMONELLA TYPHIMURIUM;                                                                        .\n',
#[Out]#  'AUTHOR    K.DIEDERICHS,W.WELTE                                                                                                 .\n',
#[Out]#  ' 1239  3  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .\n',
#[Out]#  ' 44132.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n',
#[Out]#  '  910 73.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .\n',
#[Out]#  '    6  0.5   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n']
differences = [(a, b) for a, b in zip(home, database) if a != b]
len(differences)
#[Out]# 738
len(a)
#[Out]# 137
len(zip(home, database))
#[Out]# 1269
zip(home, database)[0]
#[Out]# ('==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000 ==== DATE=15-FEB-2013     .\n',
#[Out]#  '==== Secondary Structure Definition by the program DSSP, CMBI version by M.L. Hekkelman/2010-10-21 ==== DATE=2012-03-20        .\n')
zip(home, database)[1]
#[Out]# ('REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .\n',
#[Out]#  'REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .\n')
zip(home, database)[2]
#[Out]# ('HEADER    OUTER MEMBRANE PROTEIN                  07-DEC-97   1A0S                                                             .\n',
#[Out]#  'HEADER    OUTER MEMBRANE PROTEIN                  07-DEC-97   1A0S                                                             .\n')
len(home)
#[Out]# 1269
len(differences)
#[Out]# 738
differences[0]
#[Out]# ('==== Secondary Structure Definition by the program DSSP, updated CMBI version by ElmK / April 1,2000 ==== DATE=15-FEB-2013     .\n',
#[Out]#  '==== Secondary Structure Definition by the program DSSP, CMBI version by M.L. Hekkelman/2010-10-21 ==== DATE=2012-03-20        .\n')
differences[1]
#[Out]# (' 44143.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n',
#[Out]#  ' 44132.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n')
differences[2]
#[Out]# ('  630 50.8   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n',
#[Out]#  '  627 50.6   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n')
differences[3]
#[Out]# ('    1   71 P S              0   0  104      0, 0.0    57,-0.1     0, 0.0     2,-0.1   0.000 360.0 360.0 360.0-130.6  -45.8    0.7    8.8\n',
#[Out]#  '    1   71 P S              0   0  106      0, 0.0    57,-0.1     0, 0.0     2,-0.1   0.000 360.0 360.0 360.0-130.6  -45.8    0.7    8.8\n')
differences[4]
#[Out]# ('    2   72 P G        +     0   0   29      1,-0.2    54,-3.0    55,-0.1     2,-0.3  -0.097 360.0   5.3-142.3-116.6  -43.7    1.8    5.8\n',
#[Out]#  '    2   72 P G        +     0   0   27      1,-0.2    54,-3.0    55,-0.1     2,-0.3  -0.097 360.0   5.3-142.3-116.6  -43.7    1.8    5.8\n')
differences[5]
#[Out]# ('    3   73 P F  E     -A   55   0A  51     52,-0.2     2,-0.4    -2,-0.1    -1,-0.2  -0.602  48.5-175.1 -90.0 144.8  -41.9    0.3    2.9\n',
#[Out]#  '    3   73 P F  E     -A   55   0A  50     52,-0.2     2,-0.4    -2,-0.1    -1,-0.2  -0.602  48.5-175.1 -90.0 144.8  -41.9    0.3    2.9\n')
differences[6]
#[Out]# ('    5   75 P F  E     +A   53   0A  12     -2,-0.4     2,-0.2    48,-0.2    48,-0.2  -0.953  23.9 161.8-121.9 110.0  -40.1   -5.1   -1.0\n',
#[Out]#  '    5   75 P F  E     +A   53   0A  11     -2,-0.4     2,-0.2    48,-0.2    48,-0.2  -0.953  23.9 161.8-121.9 110.0  -40.1   -5.1   -1.0\n')
exit()
