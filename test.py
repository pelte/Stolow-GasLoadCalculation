import gasLoadLibrary as gL;
import matplotlib.pyplot as plt;
import numpy as np;
import csv;

# print(gL.getWFactor(1.9,3.39))
# # print(gL.getWFactor(1,11))
# # print(gL.getWFactor(5000,11))
# # print(gL.getWFactor(5000,4.5))
# # print(gL.getWFactor(1000,200))

# print(gL.molecularSpeed(293))

# gL.displayWPlot(1.9, 3.39)

if True: print("hello")


# loading pump curves...
# GV80 = [[],[]]  # first elem is P array, second elem is S array
# with open('HeliumCellGasLoad/GV80_EH500_60Hz.csv') as pump_curve:
#     reader = csv.reader(pump_curve);
#     next(reader);       # advance pointer one row
#     for elem in reader:
#         GV80[0].append(elem[0]);
#         GV80[1].append(elem[1]);
#     #endfor
# #endwith

# fig, ax = plt.subplots(figsize = (10,4));
# ax.plot(GV80[0],GV80[1], label = "GV80 EH500 pump curve")

d = 2.5;
w = gL.beamRadius(d);
print("beam radius at a distance of {:.1f} mm from focus is: {:.2f} um".format(d,w))

# P, R, M= gL.heliumCellSimulation(100,10000, num=10000);
P,R, M = gL.heliumCellSimulation(100,1010, num = 100000);
# gL.plotCellSim(P[0], P[1], P[2], P[3])
gL.plotCellSim(P, R, M)


# fig, ax = plt.subplots(figsize = (10,4))
# ax.set_yscale("log")
# ax.set_xscale("log")


# l1, = ax.plot(P[0], P[1], marker = "o", label = "Inner Jacket");
# l2, = ax.plot(P[0], P[2], marker = "o", label = "Outer Jacket");
# l3, = ax.plot(P[0], P[3], marker = "o", label = "Lesker Chamber");
# l1, = ax.plot(P[0], P[1], label = "Inner Jacket");
# l2, = ax.plot(P[0], P[2], label = "Outer Jacket");
# l3, = ax.plot(P[0], P[3], label = "Lesker Chamber");
# ax.set_ylabel("Pressure [mbar]");
# ax.set_xlabel("Helium cell Pressure [mbar]");
# ax.set_title("3 Stage Helium Cell Pressures");

# ax.legend(handles=[l1, l2,l3])

# P = gL.heliumCellSimulation(1,10000, num=100, cell_length=10);
# gL.plotCellSim(P[0], P[1], P[2], P[3])
# fig, ax = plt.subplots(figsize = (10,4))
# ax.set_yscale("log")
# # ax.set_xscale("log")

# ax.set_ylabel("Pressure [mbar]");
# ax.set_xlabel("Helium cell Pressure [mbar]");
# ax.set_title("3 Stage Helium Cell Pressures (10 mm cell length)");
# l1, = ax.plot(P[0], P[1], label = "Inner Jacket");
# l2, = ax.plot(P[0], P[2], label = "Outer Jacket");
# l3, = ax.plot(P[0], P[3], label = "Lesker Chamber");

# ax.legend(handles=[l1, l2,l3])
# plt.show(block = True);

# testing conductance functions...
P_test = np.linspace(1E-4,1, num = 15);
C_mf = [];
C_lf = [];
T = 293;
d = .1;

for p in P_test:
	C_mf.append(gL.conductance_MF(d,T));
	C_lf.append(gL.conductance_LF(d, 0, p*2));

fig1, ax1 = plt.subplots(figsize = (10,4));
# ax1.set_yscale("log")

l4, = ax1.plot(P_test, C_mf);
l5, = ax1.plot(P_test, C_lf);

print(C_mf)

plt.show(block = True);