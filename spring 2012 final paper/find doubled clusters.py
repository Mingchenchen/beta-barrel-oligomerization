import csv

clusterlist = list()

with open("daniel's clustering info.csv", 'r') as f:
    for line in csv.reader(f):
        clusterlist += filter(lambda x: x!='', line)


cluster_count = dict()
for cluster in clusterlist:
    if cluster in cluster_count.keys():
        cluster_count[cluster] += 1
    else:
        cluster_count[cluster] = 1

for cluster, count in sorted(cluster_count.items(),
                             key=lambda x: x[1]):
    print(cluster + ': ' + str(count)) 
