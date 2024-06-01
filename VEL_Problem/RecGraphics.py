import matplotlib.pyplot as plt
import numpy as np

t, rec1, rec2, rec3 = [], [], [], []

with open("ReceiversResults.txt") as file:
    for line in file:
        tt, rec1t, rec2t, rec3t = line.split()
        t.append(float(tt))
        rec1.append(float(rec1t))
        rec2.append(float(rec2t))
        rec3.append(float(rec3t))

#with open("kek.txt") as file:
#    for line in file:
#        u.append(float(line))


# plot

plt.plot(t, rec1); name1 = 'Receiver1'
plt.plot(t, rec2); name2 = 'Receiver2'
plt.plot(t, rec3); name3 = 'Receiver3'
plt.title('Электрическая напряженность')
#plt.axis([1e-6, 1e-2, 1e-10, 1e-1])
plt.xlabel('t, с')
plt.xscale('log')

plt.ylabel('E, В/м', rotation = 0)
plt.yscale('log')
plt.grid()
#plt.xticks(rotation=-90)
plt.legend([name1, name2, name3])
plt.show()
