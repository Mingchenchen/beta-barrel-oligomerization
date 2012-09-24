# Make the modules I wrote available
import sys
sys.path.append(r'/home/davis/work stuff/beta-barrel-oligomerization/modules')

# Run the script
import os
os.chdir(r"/home/davis/work stuff/beta-barrel-oligomerization/seq id matrix mapping")
execfile(r"/home/davis/work stuff/beta-barrel-oligomerization/seq id matrix mapping/code.py")

