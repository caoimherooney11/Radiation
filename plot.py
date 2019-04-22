from firedrake import *
import matplotlib.pyplot as plt
import csv

#with open("Output/effective_conductivity.csv", 'r') as f:
with open("Output/effective_conductivity.csv", 'r') as f:
    reader = csv.reader(f, delimiter=",")
    lists = [ [] for _ in range(5) ]  
    for row in reader:
        for i in range(5):
            lists[i].append(float(row[i]))

#plt.plot(lists[0], lists[1], 'ro', lists[0], lists[2], 'bs', lists[0], lists[3], 'g^', lists[0], lists[4], 'yo')
plt.plot(lists[0], lists[1], 'ro', lists[0], lists[4], 'yo')
plt.show()

