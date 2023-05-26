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

# d = 2.5;
# w = gL.beamRadius(d);
# print("beam radius at a distance of {:.1f} mm from focus is: {:.2f} um".format(d,w))

# # P, R, M= gL.heliumCellSimulation(100,10000, num=10000);
# P,R,M = gL.heliumCellSimulation(100,10000, num = 10000);
# gL.plotCellSim(P, R, M)


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
# P_test = np.linspace(1E-4,1, num = 15);
# C_mf = [];
# C_lf = [];
# T = 293;
# d = .1;

# for p in P_test:
# 	C_mf.append(gL.conductance_MF(d,T));
# 	C_lf.append(gL.conductance_LF(d, 0, p*2));

# fig1, ax1 = plt.subplots (figsize = (10,4));
# # ax1.set_yscale("log")

# l4, = ax1.plot(P_test, C_mf);
# l5, = ax1.plot(P_test, C_lf);

# print(C_mf)

# plt.show(block = True);

# viscous flow conductance comparison
C_1 = (1E3)*gL.PI*(10**4)*(.05 + .1)/(256*(gL.MU_0*10)*1000);
# C_2 = (1E3)*gL.PI*(10**4)*(.05 + .1)*100/(256*gL.MU_0*1000);

print(C_1/1000, gL.conductance(10, .05, .1, L = 1000), gL.conductance_LF(10,.05,.1, L = 1000))

# a1, a2 = gL.conductance_MF(10,295);
# print(a1/1000, a2/1000)

M = gL.ATOMIC_M;

# molecular flow conductance comparison

C_3 = [];
C_4 = [];
LR_ratio = [];

# 2L/D > 100
print("L/a > 100");
D = 10;
L = 1000;
T = 293;

C_3.append((gL.averageSpeed(T)*100)*gL.PI*(D**3)/(12*1000*L));
C_4.append(30.48*(gL.m.sqrt(T/M))*((D/2)**3)/L);
LR_ratio.append(L/(D/2));

print(C_3[0], C_4[0]);

# 2L/D > 1.5
print("L/a > 1.5");
D = 5;
L = 500;

F_O = 11.428 *((D/2)**2)*gL.m.sqrt(T/M);
F_T = F_O/(1 + 3*L/(8*D/2));

C_3.append((gL.averageSpeed(T)*100)*gL.PI*(D**3)/(12*1000*L));
C_4.append((F_T*F_O)/(F_O + F_T));
LR_ratio.append(L/(D/2));

print(C_3[1], C_4[1]);

D = 5;
L = 10;

F_O = 11.428 *((D/2)**2)*gL.m.sqrt(T/M);
F_T = F_O/(1 + 3*L/(8*D/2));

C_3.append((gL.averageSpeed(T)*100)*gL.PI*(D**3)/(12*1000*L));
C_4.append((F_T*F_O)/(F_O + F_T));
LR_ratio.append(L/(D/2));

print(C_3[2], C_4[2]);

# 2L/D <= 1.5
print("L/a <= 1.5");
D = 5;
L = 3.75;

F_O = 11.428 *((D/2)**2)*gL.m.sqrt(T/M);
F_T = F_O/(1 + .5*(L/(D/2)));

C_3.append((gL.averageSpeed(T)*100)*gL.PI*(D**3)/(12*1000*L));
C_4.append((F_T*F_O)/(F_O + F_T));
LR_ratio.append(L/(D/2));

print(C_3[3], C_4[3]);

# fig2, ax2 = plt.subplots(figsize = (10,4));
# l1, = ax2.plot(LR_ratio, C_3);
# l2, = ax2.plot(LR_ratio, C_4);

# plt.show(block = True);

D = 5;
L = 5;

C_5 = [];
C_6 = [];
LR_ratio2 = [];

for i in range(100):
	lr_ratio = L/(D/2)
	LR_ratio2.append(lr_ratio);
	C_5.append((gL.averageSpeed(T)*100)*gL.PI*(D**3)/(12*1000*L));

	F_O = 11.428 *((D/2)**2)*gL.m.sqrt(T/M);
	if (lr_ratio >= 100):
		C_6.append(30.48*(gL.m.sqrt(T/M))*((D/2)**3)/L);
	elif lr_ratio > 1.5:
		F_T = F_O/(1 + 3*L/(8*D/2));
		C_6.append((F_T*F_O)/(F_O + F_T));
	else:
		F_T = F_O/(1 + .5*(L/(D/2))); 
		C_6.append((F_T*F_O)/(F_O + F_T));
	L = L + 2
#endfor

fig3, ax3 = plt.subplots(figsize = (10,4));
l3, = ax3.plot(LR_ratio2, C_5, label = "pfieffer eq.");
l4, = ax3.plot(LR_ratio2, C_6, label = "adjust vacuum technique eq.");


# note from the above plot that the adjusted

# testing transition formula
P_avg = (.1)/2 	#[mbar]

rho_1_sqrt	= gL.m.sqrt(M/(83.144E6 * T));
c1 		= 2*(D/2)*rho_1_sqrt/(gL.MU_0*10);
c2 		= 2.76*(D/2)*rho_1_sqrt/(gL.MU_0*10);

L = 5;
C_7 = [];
for i in range(100):
	alpha_1 	= .3926*((D/2)**4)/(gL.MU_0*10*L)
	alpha_2 	= 30476*((D/2)**3)*gL.m.sqrt(T/M)/L
	C_7.append((alpha_1*P_avg*1E3 + alpha_2*(1 + c1*P_avg*1E3)/(1 + c2*P_avg*1E3))/1000)
	L = L + 2;

l5, = ax3.plot(LR_ratio2, C_7, label = "unified equation");
ax3.legend(handles = [l3,l4,l5])

# plt.show(block = True);

# viscous conductance formula checking
P1 = .005;	# pump inlet pressure (pipe outlet)
P2 = 1;		# pipe inlet pressure

D = 5;		# 5 cm diameter
L = 50;		# .5 m initial length

C_8 		= [];	# conductance from laminar flow equation
C_9 		= [];	# conductance from unified transition eq.
LR_ratio3 	= [];	# length over radius ratio
for i in range(100):
	P_avg = (P1 + P2)/2
	LR_ratio3.append(L/(D/2));
	alpha_1 	= .3926*((D/2)**4)/(gL.MU_0*10*L)
	alpha_2 	= 30476*((D/2)**3)*gL.m.sqrt(T/M)/L
	C_9.append((alpha_1*P_avg*1E3 + alpha_2*(1 + c1*P_avg*1E3)/(1 + c2*P_avg*1E3))/1000)
	C_8.append(gL.conductance_LF(D, P1, P2, L = L));
	L = L + 2;

fig4, ax4 = plt.subplots(figsize = (10,4));
l6, = ax4.plot(LR_ratio3, C_8, label = "pfieffer eq.");
l7, = ax4.plot(LR_ratio3, C_9, label = "unified equation");
ax4.legend(handles = [l6,l7])

# checking unified transition equation over wide pressure range
dP = .005; 	# [mbar] delta pressure
P_max = 1 # [mbar] 

D = 10;
L = 1000;

P_trace = ([],[],[],[],[]);

for p in np.arange(1E-4, P_max, dP):
	# print(p)
	P_trace[0].append(p);
	P_trace[1].append(gL.conductance(D, p - dP, p, L = L));
	P_trace[2].append(gL.conductance_LF(D, p - dP, p, L = L));
	P_trace[3].append(gL.conductance_MF(D, L = L));
	P_trace[4].append(gL.Ctube_He(D, L, (p-dP/2))*.75);

fig7, ax7 = plt.subplots(figsize = (10,4))
l8, = ax7.plot(P_trace[0], P_trace[1], label = "unified conductance equation");
l9, = ax7.plot(P_trace[0], P_trace[2], label = "viscous conductance formula");
l10, = ax7.plot(P_trace[0], P_trace[3], label = "molecular conductance formula");
l11, = ax7.plot(P_trace[0], P_trace[4], label = "Leybold Formula");

ax7.legend();
ax7.set_xscale("log");

# Q = np.divide(P_trace[1],P_trace[2]);
# H = np.divide(P_trace[1],P_trace[3]);

# fig9,ax8 = plt.subplots(figsize = (10,4));
# ax8.plot(P_trace[0], Q);
# ax8.plot(P_trace[0], H);

# ax7.plot(P_trace[0], Q);

# l11, = ax7.plot(P_trace[0], Q);
# l12, = ax7.plot(P_trace[0], H);



#checking helium cell sim...
P,R,M,C = gL.heliumCellSimulation(100,10000, num = 10000);
fig5, ax5 = gL.plotCellSim(P, R, M, C)


P,R,M,C = gL.heliumCellSimulationV2(100,10000, num = 10000, cell_length = 5);
fig6, ax6 = gL.plotCellSim(P, R, M, C)

# fig5, ax5 = plt.subplots(figsize = (10,4));


plt.show(block = True);