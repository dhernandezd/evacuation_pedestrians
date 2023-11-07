import numpy as np

file_names=['../contact_force.txt', '../social_force.txt']
contact_fr = []
social_fr = []
N = 300
file_i = file_names[0]
with open(file_i, 'r') as data_file:
    for j in range(10000):
        line = data_file.readline()
        num_forces = int(line.split()[1])
        contact_fr.append(1.0*num_forces/N)
        for i in range(num_forces):
            line = data_file.readline()

file_i = file_names[1]
with open(file_i, 'r') as data_file:
    for j in range(10000):
        line = data_file.readline()
        num_forces = int(line.split()[1])
        social_fr.append(1.0*num_forces/N)
        for i in range(num_forces):
            line = data_file.readline()

np.savetxt("n_contacts.dat", np.column_stack((contact_fr, social_fr)))
