
import random
import numpy as np
import statistics
import matplotlib.pyplot as plt
import time
import datetime



start_time = time.time()

def list_txt_gen(filename, filepath, list_txt):
    with open(filepath + filename, 'w') as f:
        for i in list_txt:
            f.write("%s\n" % i)

def size_cluster(cluster, length): #calcu-lates the size of the cluster
    return float(len(cluster) / (length * length))


def find_max(object):
    max = 0
    for counter in range(len(object)):
        if max < abs(object[counter]):
            max = abs(object[counter])
    return max

def mid_point(list):
    max = 0
    min = 0
    mp = 0
    for counter in range(len(list)):
        if min > abs(list[counter]):
            min = abs(list[counter])

        if max < abs(list[counter]):
            max = abs(list[counter])

    mp_value = 0.5 * (max - min)
    for secondcounter in range(len(list)):
        if abs(list[secondcounter]) <= mp_value:
            mp = secondcounter
            break
    return mp



def En_of_S(jay, state, length, bond_vert, bond_hor):
    energy = 0
    delE = -0.005 * jay
    number = length * length
    for i in range(number):
        energy += delE * bond_vert[i] * state[i] * state[(i - length) % number]
        energy += delE * bond_hor[i] * state[i] * state[(i // length) * length + (i - 1) % length]

    return energy

def tot_mag_of_S(state, number, energy, temperature):
    tot_spin = (sum(state) / number) * np.exp((-1 * energy) / temperature)

    return tot_spin

# Set-up values :

L = 128
N = L * L

J = 1.0 # strength of mag field
T = 2.3 # temperature
del_T = 1.7 # temp range
p = 1.0 - np.exp(-1*J / T) # probability of wolff
p2 = 0 #prob of bond frustration

nsteps = 10000
ntsteps = 17
T_inc = del_T / ntsteps



mag = 0

plotT = []
avmag = []





l_row = -1
l_coloumn = -1


nbr = {i: ((i // L) * L + (i + 1) % L, (i + L) % N, (i // L) * L + (i - 1) % L, (i - L) % N) for i in range(N)} #defines a neighbourhood

tau_row = [1 for k1 in range(N)]
tau_coloumn = [1 for k2 in range(N)]

# with probability p2 assigns frustrated bonds
for x in range(N):
    if random.uniform(0.0, 1.0) < p2:
        tau_coloumn[x] = -1
    if random.uniform(0.0, 1.0) < p2:
        tau_row[x] = -1


tot_nrun = 1

step_time = time.time()


for temp in range(ntsteps):

    T += T_inc
    p = 1.0 - np.exp(-2 * J / T)
    print("duration of last cycle is : ", time.time() - step_time)
    print("temp is ", T, " while step number is ", temp, "prob is ", p)
    step_time = time.time()

    print("----------------------------------------------")

    Z = 0
    mag = 0
    S = [random.choice([1, -1]) for k in range(N)]
    for step in range(nsteps):

        E = En_of_S(J, S, L, tau_coloumn, tau_row)
        Z += np.exp((-1 * E) / T)
        # mag[step] = tot_mag_of_S(S, N, E, T)
        mag += abs(tot_mag_of_S(S, N, E, T))
        k = random.randint(0, N - 1)
        Pocket, Cluster = [k], [k]
        while Pocket != []:
            j = random.choice(Pocket)
            for nn in nbr[j]:
                if nn not in Cluster:

                    if random.uniform(0.0, 1.0) < p:
                        l_row = -10
                        l_coloumn = -5
                        if nn == (j - L) % N:
                            l_coloumn = j

                            if tau_coloumn[l_coloumn] == -1 and S[nn] == -S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                            if tau_coloumn[l_coloumn] == 1 and S[nn] == S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                        if nn == (j // L) * L + (j - 1) % L:
                            l_row = j

                            if tau_row[l_row] == 1 and S[nn] == S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                            if tau_row[l_row] == -1 and S[nn] == -S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                        if nn == (j // L) * L + (j + 1) % L:
                            l_row = nn

                            if tau_row[l_row] == 1 and S[nn] == S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                            if tau_row[l_row] == -1 and S[nn] == -S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                        if nn == (j + L) % N:
                            l_coloumn = nn

                            if tau_coloumn[l_coloumn] == -1 and S[nn] == -S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

                            if tau_coloumn[l_coloumn] == 1 and S[nn] == S[j] and nn not in Cluster:
                                # print("first", Pocket)
                                Pocket.append(nn)
                                Cluster.append(nn)
                                # print("second", Pocket)

            Pocket.remove(j)

        for j in Cluster:
            S[j] *= -1

    avmag.append(mag / Z)
    plotT.append(T)

print("----------------------------------------------------------------")


filepath = "/home/ario/PycharmProjects/Q1-git/Q1 Texts/"
filename1 = "2.3-4.0_avmag-128- 100000.txt"
filename3 = "2.3-4.0_T-128- 100000.txt"

list_txt_gen(filename1, filepath, avmag)
list_txt_gen(filename3,filepath, plotT)

end_time = time.time()

print("total runtime was : ", end_time - start_time)

end_time = time.time()
current_date_time = datetime.datetime.now()

print("total runtime was : ", end_time - start_time)

tot_time_hours = (end_time - start_time) // 3600
tot_time_mins = (end_time - start_time) / 3600 - (end_time - start_time) // 3600

#fig, ax = plt.subplots()
#plt.title('phase transition of Magnetization for L =' + str(L) + 'and steps #' + str(nsteps))
#plt.xlabel('temperature')
#plt.ylabel('absolute value of average magnetization')
#plt.plot(plotT, avmag, color='g', label='Magnetization')
#plt.axvline(x=2.269, color='red', linestyle='--')
#plt.legend()
#plt.show()

file_title = "timelog.txt"
with open('/home/ario/PycharmProjects/Q1-git/Q1 Texts/' + file_title, 'w') as f:
    f.write("%s\n" % "date and time of run" + str(current_date_time))
    f.write("%s\n" % "L =" + str(L) + "nsteps =" + str(nsteps) + "T =" + str(T))
    f.write("%s\n" % "total time is" + str(tot_time_hours) + "h &" + str(tot_time_mins) + "mins")
    f.write("%s\n" % "---------------------------------------------------------")