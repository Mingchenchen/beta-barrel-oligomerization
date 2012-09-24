# Run the script
import os
os.chdir(r"C:\cygwin\home\alex\beta-barrel-oligomerization\seq id matrix mapping")
execfile(r"C:\cygwin\home\alex\beta-barrel-oligomerization\seq id matrix mapping\code.py")

p1 = to_mat(pam1, published_bbtm_ordering)
p250 = np.linalg.matrix_power(p1, 250)
scores = np.empty(p250.shape, dtype='float')
import itertools
import math
for i, j in itertools.product(range(p250.shape[0]), range(p250.shape[1])):
	scores[i,j] = (1/0.231049) * math.log(p250[i,j] \
                                / pi_ver[published_bbtm_ordering[j]])

print(scores[0,0])
