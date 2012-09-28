import os
import sys

sys.path.append("/home/davis/work stuff/beta-barrel-oligomerization/modules")

os.chdir("/home/davis/work stuff/beta-barrel-oligomerization/zdiff module")
execfile("zdiff.py")

zdiff_align('gonnet', 'gonnet zdiff alignments')
