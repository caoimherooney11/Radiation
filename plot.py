from firedrake import *
import matplotlib.pyplot as plt
import csv
from interpolate_keff import keff


with open("Output/BigDatasets/effective_conductivity.csv", 'r') as f:
    reader = csv.reader(f, delimiter=",")
    lists = [ [] for _ in range(5) ]  
    for row in reader:
        for i in range(5):
            lists[i].append(float(row[i]))

#plt.plot(lists[0], lists[1], 'ro', lists[0], lists[2], 'bs', lists[0], lists[3], 'g^', lists[0], lists[4], 'yo')
plt.plot(lists[0], lists[1], 'ro', lists[0], keff(lists[0], lists, 1), 'y')
plt.show()
print("single value = ", keff(0.5, lists, 1))
print("matrix = ", keff(0.5, lists))
print("derivative = ", keff(0.5, lists, entry=1, derivative=1)) 
print("matrix derivative = ", keff(0.5, lists, derivative=1)) 
