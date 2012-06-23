import glob

pathnames = glob.glob('clusters/cluster*.clu')
print(pathnames)

found1 = False
found2 = False

count = 0

for pathname in pathnames:
    found1 = False
    found2 = False
    with open(pathname, 'r') as f:
        for line in f:
            if found1 and found2:
                break
            if '22.4.5' in line:
                found1 = True
            if '22.2.4' in line:
                found2 = True

    if found1 and found2:
        print(pathname)

    count +=1
    if count%10 == 0:
        print(str(count) + ' searched')


print('done')
