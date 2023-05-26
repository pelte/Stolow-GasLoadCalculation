''' Gas Load Library
by: Patrick Elten <pelten@uottawa.ca>
		started: 2023-02-13
		last updated:2023-05-26

This file will hold the functions facilitating the calculation of the gas load
target pressures.

initally created to validate the helium cell design of the HHG beamline.

bibliography:
[1] Scientific Foundations of Vacuum Technique 1st ed. 1949
[2] Building Scientific Apparatus 4th ed. 2009
[3] Rarefied gas flow through short tubes into vacuum [DOI: 10.1116/1.2830639] 2008
'''

import csv;
import math as m;
import numpy as np;
import matplotlib.pyplot as plt;
import matplotlib.ticker

# constants
# NOTE: the atom specific constants below are tuned for Helium Gas He
ATOMIC_M	= 4.003					# Atomic Mass of Helium Gas [amu]
ATOMIC_M_KG = ATOMIC_M*1.66054E-27	# Atomic Mass of Helium Gas [kg]
MOLAR_M		= .004003				# Molar Mass of He [kg/mol]
MU_0		= 0.0000197				# Dynamic Viscosity of He [Pa*s] (assumed constant)
ATOMIC_R	= 140					# Atomic Radius of He_2 [pm]
COLL_A		= .21 					# Collision Cross Section [nm^2]
# VISCOSITY	= 1.825E-5				# Dynamic viscosity of air at 20 degrees C [Pa*s]


PI			= 3.1416
AVO_NUM 	= 6.0221367E23			# Avogadro's Number
K_CONST 	= 1.38064852E-23		# Boltzmann Constant [m^2kgs^-2K^-1]
R_CONST		= 83.1446261815324		# Universal Gas Constant [L*mbar*K^-1*mol^-1]

def getWFactor(delta, in_r):
	'''
		given rarefraction parameter delta, duct length L, and aperture radius R
		return the dimensionless W factor. delta, L, R should be > 0.
		see Table 1 of "Rarefied gas flow through short tubes into vacuum" [DOI: 10.1116/1.2830639]

	This function DOES NOT EXTRAPOLATE IT ONLY INTERPOLATES. Given values outside the table range 
	will return a W factor at the nearest extreme.

		Parameters:
			delta 	(float): rarefraction parameter of gas
			in_r 	(float): length over radius ratio

		Returns:
			W 		(float): W factor based on Table 1 in [3]
	'''
	ratioList = [0, 0.1, 0.5, 1, 5, 10];	# elements in this list are columns of Table 1
	deltaList = [0, .1, .5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000];	# elements in this list are rows of Table 1

	# determine the column number for the L/r value closest to but LOWER than the given L/R ratio
	if in_r >= ratioList[-1]: col = len(ratioList) - 1;
	else:
		col = 0;
		for itr in enumerate(ratioList):
			if in_r <= itr[1]:
				col = itr[0] - 1;
				break;

	# determine the row number for the delta value closest to but LOWER than the given delta value
	if delta >= deltaList[-1]: row = len(deltaList) - 1;
	else:
		row = 0;
		for itr in enumerate(deltaList):
			if delta <= itr[1]:
				row = itr[0] - 1;
				break;

	# print("(row, column) = ({:d}, {:d})".format(row, col))

	# this list will hold the W values read from the csv, used to interpolate
	W = [0,0,0,0]

	with open('resources/Table1.csv') as table1:
		reader = csv.reader(table1);
		i = 0;
		while(i < row):
			next(reader);
			i = i +1;

		current = next(reader)
		W[0] = float(current[col]);
		if col == len(ratioList) - 1:		# if at extreme merge W_1 and W_0
			W[1] = W[0];
		else:
			W[1] = float(current[col + 1]);

		if i != len(deltaList) -1: current = next(reader);
		W[2] = float(current[col]);
		if col == len(ratioList) - 1:		# if at extreme, merge W_3 and W_2
			W[3] = W[2];
		else:
			W[3] = float(current[col + 1]);

	# print("W_1: {:.3f}\t W_2: {:.3f}\t W_3: {:.3f}\t W_4: {:.3f}".format(W[0],W[1],W[2],W[3]))
	if col == len(ratioList) - 1:
		ratio_factor = 0
	else:
		ratio_factor = (in_r - ratioList[col])/(ratioList[col+1] - ratioList[col])
	# print(ratio_factor)

	W_prime = [0,0];
	W_prime[0] = W[0] + (W[1] - W[0])*ratio_factor;
	W_prime[1] = W[2] + (W[3] - W[2])*ratio_factor;

	# print(W_prime[0]);
	# print(W_prime[1]);

	if row == len(deltaList) - 1:
		delta_factor = 0
	else:
		delta_factor = (delta - deltaList[row])/(deltaList[row + 1] - deltaList[row])

	retval = (W_prime[0] + (W_prime[1] - W_prime[0])*delta_factor)

	return retval;
#endfunction

def displayWPlot(delta, in_r):
	'''
		given rarefraction parameter delta, duct length L, and aperture radius R
		computes and plots the W factor. see function getWFactor(delta, in_r)

			Parameters:
				delta 	(float): rarefraction parameter for gas
				in_r 	(float): length over radius ratio of tube
			Returns:
				(fig, ax)
				fig 	(figure): figure reference for W Factor plot
				ax 		(Axes): Axes reference for W Factor plot
	'''
	ratioList = [0, 0.1, 0.5, 1, 5, 10];	# elements in this list are columns of Table 1
	deltaList = [0, .1, .5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000];	# elements in this list are rows of Table 1

	lines = [[],[],[],[],[],[]];
	# print("\n")

	with open('resources/Table1.csv') as table1:
		reader = csv.reader(table1);
		while reader.line_num < 14:
			current = next(reader);
			for i in range(6):
				elem = str(current[i]).strip()
				lines[i].append(float(elem))
				# print(current[i])

	fig, ax = plt.subplots(figsize = (10,4))
	ax.set_xscale("log")
	ax.set_xticks(deltaList)
	ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
	ax.set_ylabel("W factor")
	ax.set_xlabel("Rarefaction Parameter [$\delta$]")
	ax.set_title("W factor with lines of constant L/R [1]");

	l1, = ax.plot(deltaList, lines[0], label = "L/R = 0");
	l2, = ax.plot(deltaList, lines[1], label = "L/R = .1")
	l3, = ax.plot(deltaList, lines[2], label = "L/R = .5")
	l4, = ax.plot(deltaList, lines[3], label = "L/R = 1")
	l5, = ax.plot(deltaList, lines[4], label = "L/R = 5")
	l6, = ax.plot(deltaList, lines[5], label = "L/R = 10")

	ax.legend(handles=[l1,l2,l3,l4,l5,l6])

	retval = getWFactor(delta, in_r);

	ax.plot(delta, retval, color = 'black', marker = 'x', markersize=12, 
	label = "({:.2f},{:1f})".format(delta,retval));
	ax.annotate("({:.2f},{:.3f})".format(delta,retval),(delta,retval),
						textcoords = "offset points", xytext = (-30,-20))
	plt.show(block = True)

	return(fig,ax)
#endfunction

def massflow_lim(R, P_0, v_0):
	'''
	function returns the mass flow through an orifice at the free-molecular limit (rarefraction parameter = 0)
	Citation: [3, pg 229, eq. (3)]

		Parameters:
			R 	():
			P_0 ():
			v_0 ():

		Returns:
			M_0' (float): mass flow through an orifice at the free-molecular limit [kg/s]

	'''
	return m.sqrt(PI)*R*R*P_0/v_0;
#endfunction

def molecularSpeed(T = 293):
	'''given gas temperature returns the most probable molecular speed
		Citation [3, pg. 229]

		Parameters:
			T 	(int): gas temperature, default = 293 [K]
		Returns:
			v_0 (float): most probable molecular speed [m/s]
	'''
	return m.sqrt(2*K_CONST*T/ATOMIC_M_KG);
#endfunction

def rarefractionParameter(R, P, v_0):
	'''returns the rarefraction Parameter of the gas through the short tube
	Citation [3, pg 229, eq. (1)]

		Parameters:
			R 		(float): aperture radius [m]
			P 		(float): reference pressure [Pa]
			v_0 	(float): molecular speed [m/s]
		Returns:
			delta 	(float): dimensionless rarefraction parameter delta
	'''
	return (R*P)/(MU_0*v_0);
#endfunction

def massFlow_to_throughPut(m, T = 293):
	'''converts mass flow rate to through put rate using ideal gas law
		Parameters:
			m 	(float): mass flow [kg/s]
			T 	(int): gas temperature, default = 293
		Returns:
			Q 	(float): gas throughput [mbar*L/s]
	'''
	return m*R_CONST*T/MOLAR_M;
#endfunction

def averageSpeed(T = 293):
	'''returns the average speed of gas at the given temperature in m/s
	Citation [2, pg 94, eq. (3.1)]

		Parameters:
			T (int): gas temperature, dfault = 293 [K]
		Returns:
			v_bar (float): average speed [m/s]
	'''
	# return m.sqrt(8*K_CONST*T/(PI*ATOMIC_M_KG));
	return m.sqrt(8*(R_CONST/10)*T/(PI*MOLAR_M));
#endfunction

def conductance_MF(D,T = 293, L = 1):
	'''
		Conductance (in L/s) of a round pipe experiencing molecular flow. Default length of 1 cm.
		For small lengths, (L/(D/2) < 100) end effects are considered.
		Citation [1, pg 100-101]

			Parameters:
				D 	(float): tube internal diameter [cm]
				T 	(int): gas temperature, default = 293 [K]
				L 	(float): tube length, default = 1 [cm]
			Returns:
				C 	(float): tube conductance [L/s]
	'''
	# return (2.6E-4*(averageSpeed(T)*100)*((D)**3));		# [2, pg. 97, eq. (3.19)]
	if L == 1: 
		return (averageSpeed(T)*100)*PI*(D**3)/(12*1000*L);		# pfeiffer site
	#endif

	# taking into account end effects (for a short tube):
	lr_ratio = L/(D/2);
	F_0 = 11.428*((D/2)**2)*m.sqrt(T/ATOMIC_M);
	if(lr_ratio >= 100):
		return (30.48*(m.sqrt(T/ATOMIC_M))*((D/2)**3)/L);
	elif lr_ratio > 1.5:
		F_T = F_O/(1 + 3*L/(8*D/2));
		return ((F_T*F_O)/(F_O + F_T));
	else:
		F_T = F_O/(1 + .5*(L/(D/2))); 
		return ((F_T*F_O)/(F_O + F_T));
	#endif
#endfunction

def conductance_LF(D, P1, P2, L = 1):
	'''
		Conductance (in L/s) of a round pipe experiencing laminar viscous flow. Default length of 1 cm.
		Citation [1, pg. 86, eq. (5)]

			Parameters:
				D 	(float): tube internal diameter [cm]
				P1 	(float): pressure value at input [mbar]
				P2 	(float): pressure value at output [mbar]
				L 	(float): tube length, default = 1 [cm]
			Returns:
				C 	(float): tube conductance [L/s]
	'''
	return PI*(D**4)*(P1 + P2)/(256*MU_0*10*L);		# citation [1, pg. 86, eq. (5)]
#endfunction

def conductance(D, P1, P2, T = 293, L = 1):
	'''
		Conductance (in L/s) of a round pipe experiencing either molecular flow, laminar viscous flow or 
		within the transition regime. Default length of 1 cm.
		Citation [1, pg. 111, eq. (1) (1.2) (2.8) (2a) (3a)]

			Parameters:
				D 	(float): tube internal diameter [cm]
				P1 	(float): pressure value at input [mbar]
				P2 	(float): pressure value at output [mbar]
				T 	(int): gas temperature, default = 293 [K]
				L 	(float): tube length, default = 1 [cm]
			Returns:
				C 	(float): tube conductance [L/s]
	'''
	P_avg = (P1 + P2)/2;
	
	rho_1_sqrt	= m.sqrt(ATOMIC_M/(83.144E6 * T));
	c1 			= 2*(D/2)*rho_1_sqrt/(MU_0*10);
	c2 			= 2.76*(D/2)*rho_1_sqrt/(MU_0*10);
	alpha_1 	= .3926*((D/2)**4)/(MU_0*10*L)
	alpha_2 	= 30476*((D/2)**3)*m.sqrt(T/ATOMIC_M)/L
	return ((alpha_1*P_avg*1E3 + alpha_2*(1 + c1*P_avg*1E3)/(1 + c2*P_avg*1E3))/1000)
#endfunction

def Ctube_He (d,l,Pav):
    '''
    Knudsen flow @20degC, helium 
    (formula for rarified air from Leybold's website scaled to He, see factors sc1 and sc2)
    Conductance for a straight pipe, which is not too short, of length l, with a circular cross section of diameter d 
    for the laminar, Knudsen and molecular flow ranges, valid for air at 20 °C (Knudsen equation):
    https://content.leybold.com/en/knowledge/what-is-the-conductance-for-piping-and-orifices
    d = Pipe inside diameter in cm
    l = Pipe length in cm (l ≥ 10 d)
    Pav=p1/2+p2/2
    p1 = Pressure at start of pipe (along the direction of flow) in mbar
    p2 = Pressure at end of pipe (along the direction of flow) in mbar
    C[l/s]
    '''
    Pav=Pav/0.75 #[torr->mbar]
    sc1 = 185/199 #based on scoles book
    sc2 = 1368e2/508e2 #=nu(He)/nu(Air) based on BSI book
    #C = 135*d**4/l*Pav + 12.1*d**3/l* (1+192*d*Pav)/(1+237*d*Pav) #<-unscaled, for air
    C = sc1 * 135*d**4/l*Pav + sc2 * 12.1*d**3/l* (1+192*d*Pav)/(1+237*d*Pav)
    return C
#endfunction

def beamRadius (z, f = 1E6, lmbd = 2, r_lens = 15E3):
	'''
		Given a distance from the focus, returns the beam radius of focused beam.
		Citation [2, pg. 180, eq. 4.127]

			Parameters:
				z 		(float): distance from focus [mm]
				f 		(float): focal length, default = 1E6 [um]
				lmbd	(float): wavelength of light, default = 2 [um]
				r_lens 	(float): radius of beam on focusing lens, default = 15E3 [um]
			Returns:
				w 		(float): radius of beam at distance z from focus [um]
	'''
	# using beam diameter at lens:
	w_0 = (f*lmbd*1)/(PI*r_lens);
	# print("w_0 {:.2f} um".format(w_0))

	# rayleigh range [um];
	z_R = (PI*(w_0**2))/lmbd;
	# print(w_0*)
	# print("z_R {:.2f} um".format(z_R))

	return w_0*m.sqrt(1+((z*1000/z_R)**2));
#endfunction

def meanFreePath(P, T = 293):
	'''	returns the mean free path of the gas at a given temperature and pressure.
	Citation [2, pg. 94, eq. (3.4)]

		Parameters:
			T 		(int): gas temperature, default = 293 [K]
			P 		(float): gas pressure [mbar]
		Returns:
			lmbd	(float): mean free path [cm]
	'''
	return ((R_CONST*T/10)/(m.sqrt(2)*(COLL_A)*(1E-18)*AVO_NUM*(P*100)))*100;
#endfunction

def heliumCellSimulation(P_min, P_max, cell_length=5, num=24):
	'''
	This summary program performs an iterative simulation for the high harmonic generation helium cell
	to determine the steady state pressures inside the inner jacket, outer jacket, and lesker
	chamber while sweeping through a range of helium cell pressures as determined
	by P_min, P_max. This program attempts to determine regime by comparing the mean free path
	of the Helium gas and a characteristic dimension for each pumped stage. The corresponding
	conductance equation is then used to determine the effective pumping rate.

		Parameters:
			P_min 		(float): minimum helium cell pressure for simulation [mbar]
			P_mix 		(float): maximum helium cell pressure for simulation [mbar]
			cell_length	(float): simulated helium cell length, default = 5 [mm]
			num 		(int): integer value for the number of steps to take between P_min and
						P_max, default = 24
		Returns:
			(P, R, M, C)

			P = (P_0, P_1, P_2, P_3)
				P_0 	(array(float)): helium cell pressures [mbar]
				P_1 	(array(float)): inner jacket pressures [mbar]
				P_2 	(array(float)): outer jacket pressures [mbar]
				P_3 	(array(float)): lesker chamber pressures [mbar]

			R = (regime_1, regime_2, regime_3)
				regime_1 (array(int)): regime trace for inner jacket(1 for viscous flow, 0 for molecular)
				regime_2 (array(int)): regime trace for outer jacket(1 for viscous flow, 0 for molecular)
				regime_3 (array(int)): regime trace for lesker chamber(1 for viscous flow, 0 for molecular)

			M = (mfp_1, mfp_2, mfp_3)
				mfp_1 (array(float)): mean free path for each helium cell pressure point for inner jacket [cm]
				mfp_2 (array(float)): mean free path for each helium cell pressure point for outer jacket [cm]
				mfp_3 (array(float)): mean free path for each helium cell pressure point for lesker chamber [cm]

			C = (F_1, F_2, F_3)
				F_1 (array(float)): effective conductance of inner jacket pumping system [L/s]
				F_2 (array(float)): effective conductance of outer jacket pumping system [L/s]
				F_3 (array(float)): effective conductance of lesker chamber pumping system [L/s]
	'''
	T = 295;				# gas temperature [K] (assumed constant)
	tape_t = .0036*25.4;	# tape thicknmess [mm]

	P_0 = np.linspace(P_min, P_max, num=num, endpoint=True);
	P_1 = [];				# inner jacket pressure [mbar]
	P_2 = [];				# outer jacket pressure [mbar]
	P_3 = [];				# lesker chamber pressure [mbar]

	F_1 = [];
	F_2 = [];
	F_3 = [];

	regime_1 = [];			# flow regime marker, 0 for molecular, 1 for viscous
	regime_2 = [];	
	regime_3 = [];

	mfp_1 = [];				# mean free path trace
	mfp_2 = [];
	mfp_3 = [];


	v_0 = molecularSpeed(T);

	# replace this with spot size calculator!
	R_0 = 2*beamRadius(cell_length/2)/ 1E3;
	L_0 = tape_t;
	LR_ratio_0 = L_0/R_0;
	S_initial_1 = 130;

	dF_1 = 10           # Diameter of pumping forelines [cm]
	lF_1 = 1000         # Length of pumping forelines [cm]

	R_1 = 2*beamRadius(.5*25.4)/ 1E3;
	L_1 = tape_t;
	LR_ratio_1 = L_1/R_1;
	S_initial_2 = 140;

	dF_2 = 10           # Diameter of pumping forelines [cm]
	lF_2 = 1000         # Length of pumping forelines [cm]

	R_2 = 3;				# conductance hole radius between outer jacket and chamber [mm]
	L_2 = 10*R_2;
	LR_ratio_2 = L_2/R_2;
	S_initial_3 = 140;

	dF_3 = 15          # Diameter of pumping forelines [cm]
	lF_3 = 1000        # Length of pumping forelines [cm]

	# loading pump curves...
	GV80 = [[],[]]  # first elem is P array, second elem is S array
	with open('resources/GV80_EH500_60Hz.csv') as pump_curve:
	    reader = csv.reader(pump_curve);
	    next(reader);       # advance pointer one row
	    for elem in reader:
	        GV80[0].append(float(elem[0]));
	        GV80[1].append(float(elem[1]));
	    #endfor
	#endwith

	EPX500 = [[],[]]  # first elem is P array, second elem is S array
	with open('resources/EPX500LE.csv') as pump_curve:
	    reader = csv.reader(pump_curve);
	    next(reader);       # advance pointer one row
	    for elem in reader:
	        EPX500[0].append(float(elem[0]));
	        EPX500[1].append(float(elem[1]));
	    #endfor
	#endwith

	i = 0;

	# The following variables hold the previous state values for each chamber
	# (pumping speed, chamber pressure, pump inlet pressure)

	# inner jacket
	delta_0 = rarefractionParameter(R_0/1E3, P_0[0]*100, v_0);
	W_0 = getWFactor(delta_0, LR_ratio_0);
	M_0 = massflow_lim(R_0/1E3, P_0[0]*100, v_0)*W_0;
	Q_0 = massFlow_to_throughPut(M_0,T);
	prev_state_1 = (S_initial_1, (2*Q_0)/S_initial_1, (2*Q_0)/S_initial_1);

	# outer jacket
	delta_1 = rarefractionParameter(R_1/1E3, prev_state_1[2]*100, v_0);
	W_1 = getWFactor(delta_1, LR_ratio_1);
	M_1 = massflow_lim(R_1/1E3, prev_state_1[2]*100, v_0)*W_1;
	Q_1 = massFlow_to_throughPut(M_1,T);
	prev_state_2 = (S_initial_2, (2*Q_1)/S_initial_2, (2*Q_1)/S_initial_2);

	# lesker chamber
	delta_2 = rarefractionParameter(R_2/1E3, prev_state_2[2]*100, v_0);
	W_2 = getWFactor(delta_2, LR_ratio_2);
	M_2 = massflow_lim(R_2/1E3, prev_state_2[2]*100, v_0)*W_2;
	Q_2 = massFlow_to_throughPut(M_2,T);
	prev_state_3 = (S_initial_3, (Q_2)/S_initial_3, (Q_2)/S_initial_3);


	for p in P_0:
		# print("\nhelium cell pressure... {:.2f} mbar".format(p));
		##############################################################################################
		# helium cell
		# print("\ncell to inner jacket...")

		delta_0 = rarefractionParameter(R_0/1E3, p*100, v_0);
		W_0 = getWFactor(delta_0, LR_ratio_0);
		M_0 = massflow_lim(R_0/1E3, p*100, v_0)*W_0;
		Q_0 = massFlow_to_throughPut(M_0,T);

		err = 0;
		first = True;
		regime_flag = 0;
		d = 1*2.54;		# viscous transition around mean free path ~= 1 inches
		low_d 		= d/3;
		high_d 		= 3*d;
		while(first or err > 2):
			if(first):
				S_1 = prev_state_1[0];
				P1 = prev_state_1[1];
				if (P1 < 1E-3):
					P1 = 1E-3
				#endif
				P_pump1 = prev_state_1[2];
				first = False;
	        #endif

	        #1. calculate conductance
			# if(meanFreePath(P1, T) > 2*2.54):	# viscous edge at mean free path > 2 inches
			# 	# print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P1, T), P1));
			# 	C_1 = conductance_MF(dF_1, T)/lF_1;	            # molecular flow
			# 	regime_flag = 0;
			# else:
			# 	# print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P1, T), P1));
			# 	C_1 = conductance_LF(dF_1, P1, P_pump1)/lF_1;   # viscous flow
			# 	regime_flag = 1;

			#1. calculate conductance
			x = meanFreePath(P1, T);
			if(x > high_d):	
				C_1 = conductance_MF(dF_1, T)/lF_1;	            # molecular flow
				# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
				# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
				# print("x... {:.2f}".format(x));
				frac = (x - low_d)/(high_d - low_d);
				C_1 = (conductance_MF(dF_1, T)/lF_1)*(frac) + (conductance_LF(dF_1, P1, P_pump1)/lF_1)*(1-frac);
				# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
					# .format(conductance_MF(dF_2, T)/lF_2, conductance_LF(dF_2, P2, P_pump2)/lF_2, C_2));
				# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
				C_1 = conductance_LF(dF_1, P1, P_pump1)/lF_1;   # viscous flow
				# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2))				
				regime_flag = 1;
	        # # endif
			

			# C_1 = conductance(dF_1, P1, P_pump1, L = lF_1);

	        #2. calculate Seff
			S_eff1 = 1/((1/S_1) + (1/C_1));

			#3. calculate chamber pressure
			P1 = (2*Q_0)/S_eff1;

	        #4. calculate pressure at pump inlet
			P_pump1 = P1 - (2*Q_0)/C_1;

	        #5. calculate error and update pumping speed
			err = abs((S_1 - np.interp(P_pump1, GV80[0], GV80[1], left = 0, right = 0))*100/S_1);
			S_1 = (S_1 + np.interp(P_pump1, GV80[0], GV80[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable
		prev_state_1 = (S_1, P1, P_pump1);

		# inner jacket pumping
		regime_1.append(regime_flag)
		mfp_1.append(meanFreePath(P1, T));
		P_1.append(P1);
		F_1.append(C_1);

		##############################################################################################
		# inner jacket to outer jacket
		# print(i)
		# print(P_1)
		# print(P_1)
		# print("\ninner jacket to outer...")

		delta_1 = rarefractionParameter(R_1/1E3, P_1[i]*100, v_0);
		W_1 = getWFactor(delta_1, LR_ratio_1);
		M_1 = massflow_lim(R_1/1E3, P_1[i]*100, v_0)*W_1;
		Q_1 = massFlow_to_throughPut(M_1,T);
		
		err = 0;
		first = True;
		regime_flag = 0;
		d = 4*2.54;		# viscous transition around mean free path ~= 4 inches
		low_d 		= d/3;
		high_d 		= 3*d;
		while(first or err > 2):
			if(first):
				S_2 = prev_state_2[0];
				P2 = prev_state_2[1];
				# if (P2 < 1E-5):
				# 	P2 = 1E-5
				#endif
				P_pump2 = prev_state_2[2];
				first = False;
	        #endif

	        # #1. calculate conductance
			# if(meanFreePath(P2, T) > 8*2.54):	
			# 	C_2 = conductance_MF(dF_2, T)/lF_2;	            # molecular flow
			# 	print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(P2, T), P2, C_2));
			# 	regime_flag = 0;
			# else:
			# 	C_2 = conductance_LF(dF_2, P2, P_pump2)/lF_2;   # viscous flow
			# 	print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(P2, T), P2, C_2))				
			# 	regime_flag = 1;
	        # # endif

	        #1. calculate conductance
			x = meanFreePath(P2, T);
			if(x > high_d):	
				C_2 = conductance_MF(dF_2, T)/lF_2;	            # molecular flow
				# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
				# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
				# print("x... {:.2f}".format(x));
				frac = (x - low_d)/(high_d - low_d);
				C_2 = (conductance_MF(dF_2, T)/lF_2)*(frac) + (conductance_LF(dF_2, P2, P_pump2)/lF_2)*(1-frac);
				# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
					# .format(conductance_MF(dF_2, T)/lF_2, conductance_LF(dF_2, P2, P_pump2)/lF_2, C_2));
				# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
				C_2 = conductance_LF(dF_2, P2, P_pump2)/lF_2;   # viscous flow
				# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2))				
				regime_flag = 1;
	        # # endif

	        #2. calculate Seff
			S_eff2 = 1/((1/S_2) + (1/C_2));

			#3. calculate chamber pressure
			P2 = (2*Q_1)/S_eff2;

	        #4. calculate pressure at pump inlet
			P_pump2 = P2 - (2*Q_1)/C_2;

	        #5. calculate error and update pumping speed
			err = abs((S_2 - np.interp(P_pump2, EPX500[0], EPX500[1], left = 0, right = 0))*100/S_2);
			S_2 = (S_2 + np.interp(P_pump2, EPX500[0], EPX500[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable...
		prev_state_2 = (S_2, P2, P_pump2);

		# outer jacket pumping
		regime_2.append(regime_flag);
		mfp_2.append(meanFreePath(P2, T));
		P_2.append(P2);
		F_2.append(C_2);

		##############################################################################################
		# outer jacket to Lesker Chamber
		# print("\nouter jacket to lesker...")
		delta_2 = rarefractionParameter(R_2/1E3, P_2[i]*100, v_0);
		W_2 = getWFactor(delta_2, LR_ratio_2);
		M_2 = massflow_lim(R_2/1E3, P_2[i]*100, v_0)*W_2;
		Q_2 = massFlow_to_throughPut(M_2,T);

		err = 0;
		first = True;
		regime_flag = 0;
		d = 16*2.54;
		low_d 		= d/3;
		high_d		= d*3;
		while(first or err > 2):
			if(first):
				S_3 = prev_state_3[0];
				P3 = prev_state_3[1];
				# if (P3 < 1E-5):
				# 	P3 = 1E-5
				#endif
				P_pump3 = prev_state_3[2];
				first = False;

	        #1. calculate conductance
			# if(meanFreePath(P3, T) > 45*2.54):		# viscous edge at mean free path > 45 inches
			# 	# print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P3, T), P3));
			# 	C_3 = conductance_MF(dF_3, T)/lF_3;	            # molecular flow
			# 	regime_flag = 0;
			# else:
			# 	# print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P3, T), P3));
			# 	C_3 = conductance_LF(dF_3, P3, P_pump3)/lF_3;   # viscous flow
			# 	regime_flag = 1;	        
			# #endif

			x = meanFreePath(P3, T);
			if(x > high_d):	
				C_3 = conductance_MF(dF_3, T)/lF_3;	            # molecular flow
				# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P3, C_3));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
				# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
				# print("x... {:.2f}".format(x));
				frac = (x - low_d)/(high_d - low_d);
				C_3 = (conductance_MF(dF_3, T)/lF_3)*(frac) + (conductance_LF(dF_3, P3, P_pump3)/lF_3)*(1-frac);
				# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
					# .format(conductance_MF(dF_3, T)/lF_3, conductance_LF(dF_3, P3, P_pump3)/lF_3, C_3));
				# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
				C_3 = conductance_LF(dF_3, P3, P_pump3)/lF_3;   # viscous flow
				# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P3, C_3))				
				regime_flag = 1;
	        # # endif

	        #2. calculate Seff
			S_eff3 = 1/((1/S_3) + (1/C_3));

			#3. calculate chamber pressure
			P3 = (Q_2)/S_eff3;

	        #4. calculate pressure at pump inlet
			P_pump3 = P3 - (Q_2)/C_3;

	        #5. calculate error and update pumping speed
			err = abs((S_3 - np.interp(P_pump3, EPX500[0], EPX500[1], left = 0, right = 0))*100/S_3);
			S_3 = (S_3 + np.interp(P_pump3, EPX500[0], EPX500[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable
		prev_state_3	 = (S_3, P3, P_pump3);

		# Lesker Chamber pumping
		regime_3.append(regime_flag);
		mfp_3.append(meanFreePath(P3, T));
		P_3.append(P3);
		F_3.append(C_3);

		i +=1;
	#endfor

	# for i in range(len(P_0)):
	# 	print("({:.1f}, {:.3e}, {:d}, {:.3e})".format(P_0[i], P_1[i], regime_1[i], mfp_1[i]));
	#endfor
	# print(P_0)
	# print(P_1);
	# print(regime_1);
	# print(mfp_1);

	return((P_0, P_1, P_2, P_3), (regime_1, regime_2, regime_3), (mfp_1, mfp_2, mfp_3), (F_1, F_2, F_3));
#endfunction

def heliumCellSimulationV2(P_min, P_max, cell_length=5, num=24):
	'''
	This summary program performs an iterative simulation for the high harmonic generation helium cell
	to determine the steady state pressures inside the inner jacket, outer jacket, and lesker
	chamber while sweeping through a range of helium cell pressures as determined
	by P_min, P_max. Gas flow regime is tracked by comparing the mean free path of the gas with the
	characteristic dimension of each pumping stage. Conductance is calculated using the function
	conductance() over all regimes.

		Parameters:
			P_min 		(float): minimum helium cell pressure for simulation [mbar]
			P_mix 		(float): maximum helium cell pressure for simulation [mbar]
			cell_length	(float): simulated helium cell length, default = 5 [mm]
			num 		(int): integer value for the number of steps to take between P_min and
						P_max, default = 24
		Returns:
			(P, R, M, C)

			P = (P_0, P_1, P_2, P_3)
				P_0 	(array(float)): helium cell pressures [mbar]
				P_1 	(array(float)): inner jacket pressures [mbar]
				P_2 	(array(float)): outer jacket pressures [mbar]
				P_3 	(array(float)): lesker chamber pressures [mbar]

			R = (regime_1, regime_2, regime_3)
				regime_1 (array(int)): regime trace for inner jacket(1 for viscous flow, 0 for molecular)
				regime_2 (array(int)): regime trace for outer jacket(1 for viscous flow, 0 for molecular)
				regime_3 (array(int)): regime trace for lesker chamber(1 for viscous flow, 0 for molecular)

			M = (mfp_1, mfp_2, mfp_3)
				mfp_1 (array(float)): mean free path for each helium cell pressure point for inner jacket [cm]
				mfp_2 (array(float)): mean free path for each helium cell pressure point for outer jacket [cm]
				mfp_3 (array(float)): mean free path for each helium cell pressure point for lesker chamber [cm]

			C = (F_1, F_2, F_3)
				F_1 (array(float)): effective conductance of inner jacket pumping system [L/s]
				F_2 (array(float)): effective conductance of outer jacket pumping system [L/s]
				F_3 (array(float)): effective conductance of lesker chamber pumping system [L/s]
	'''
	T = 295;				# gas temperature [K] (assumed constant)
	tape_t = .0036*25.4;	# tape thicknmess [mm]

	P_0 = np.linspace(P_min, P_max, num=num, endpoint=True);
	P_1 = [];				# inner jacket pressure [mbar]
	P_2 = [];				# outer jacket pressure [mbar]
	P_3 = [];				# lesker chamber pressure [mbar]

	F_1 = [];
	F_2 = [];
	F_3 = [];

	regime_1 = [];			# flow regime marker, 0 for molecular, 1 for viscous
	regime_2 = [];	
	regime_3 = [];

	mfp_1 = [];				# mean free path trace
	mfp_2 = [];
	mfp_3 = [];


	v_0 = molecularSpeed(T);

	# replace this with spot size calculator!
	R_0 = 2*beamRadius(cell_length/2)/ 1E3;
	# print(R_0)
	L_0 = tape_t;
	LR_ratio_0 = L_0/R_0;
	S_initial_1 = 130;

	# inner jacket forelines
	dF_1 = 10           # Diameter of pumping forelines [cm]
	lF_1 = 1000         # Length of pumping forelines [cm]

	R_1 = 2*beamRadius(.5*25.4)/ 1E3;
	# print(R_1)
	L_1 = tape_t;
	LR_ratio_1 = L_1/R_1;
	S_initial_2 = 140;

	# outer jacket forelines
	dF_2 = 10           # Diameter of pumping forelines [cm]
	lF_2 = 1000         # Length of pumping forelines [cm]

	R_2 = 3;				# conductance hole radius between outer jacket and chamber [mm]
	# print(R_2);
	L_2 = 10*R_2;
	LR_ratio_2 = L_2/R_2;
	S_initial_3 = 140;

	# lesker chamber forelines
	dF_3 = 15          # Diameter of pumping forelines [cm]
	lF_3 = 1000        # Length of pumping forelines [cm]

	# loading pump curves...
	GV80 = [[],[]]  # first elem is P array, second elem is S array
	with open('resources/GV80_EH500_60Hz.csv') as pump_curve:
	    reader = csv.reader(pump_curve);
	    next(reader);       # advance pointer one row
	    for elem in reader:
	        GV80[0].append(float(elem[0]));
	        GV80[1].append(float(elem[1]));
	    #endfor
	#endwith

	EPX500 = [[],[]]  # first elem is P array, second elem is S array
	with open('resources/EPX500LE.csv') as pump_curve:
	    reader = csv.reader(pump_curve);
	    next(reader);       # advance pointer one row
	    for elem in reader:
	        EPX500[0].append(float(elem[0]));
	        EPX500[1].append(float(elem[1]));
	    #endfor
	#endwith

	i = 0;

	# The following variables hold the previous state values for each chamber
	# (pumping speed, chamber pressure, pump inlet pressure)

	# inner jacket
	delta_0 = rarefractionParameter(R_0/1E3, P_0[0]*100, v_0);
	W_0 = getWFactor(delta_0, LR_ratio_0);
	M_0 = massflow_lim(R_0/1E3, P_0[0]*100, v_0)*W_0;
	Q_0 = massFlow_to_throughPut(M_0,T);
	prev_state_1 = (S_initial_1, (2*Q_0)/S_initial_1, (2*Q_0)/S_initial_1);

	# outer jacket
	delta_1 = rarefractionParameter(R_1/1E3, prev_state_1[2]*100, v_0);
	W_1 = getWFactor(delta_1, LR_ratio_1);
	M_1 = massflow_lim(R_1/1E3, prev_state_1[2]*100, v_0)*W_1;
	Q_1 = massFlow_to_throughPut(M_1,T);
	prev_state_2 = (S_initial_2, (2*Q_1)/S_initial_2, (2*Q_1)/S_initial_2);

	# lesker chamber
	delta_2 = rarefractionParameter(R_2/1E3, prev_state_2[2]*100, v_0);
	W_2 = getWFactor(delta_2, LR_ratio_2);
	M_2 = massflow_lim(R_2/1E3, prev_state_2[2]*100, v_0)*W_2;
	Q_2 = massFlow_to_throughPut(M_2,T);
	prev_state_3 = (S_initial_3, (Q_2)/S_initial_3, (Q_2)/S_initial_3);

	for p in P_0:
		# print("\nhelium cell pressure... {:.2f} mbar".format(p));
		##############################################################################################
		# helium cell
		# print("\ncell to inner jacket...")

		delta_0 = rarefractionParameter(R_0/1E3, p*100, v_0);
		W_0 = getWFactor(delta_0, LR_ratio_0);
		M_0 = massflow_lim(R_0/1E3, p*100, v_0)*W_0;
		Q_0 = massFlow_to_throughPut(M_0,T);

		err = 0;
		first = True;
		regime_flag = 0;
		d = 1*2.54;		# viscous transition around mean free path ~= 1 inches
		low_d 		= d/3;
		high_d 		= 3*d;
		while(first or err > 2):
			if(first):
				S_1 = prev_state_1[0];
				P1 = prev_state_1[1];
				if (P1 < 1E-3):
					P1 = 1E-3
				#endif
				P_pump1 = prev_state_1[2];
				first = False;
	        #endif

	        #1. calculate conductance
			x = meanFreePath(P1, T);
			if(x > high_d):	
				# C_1 = conductance_MF(dF_1, T)/lF_1;	            # molecular flow
				# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
				# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
				# print("x... {:.2f}".format(x));
				# frac = (x - low_d)/(high_d - low_d);
				# C_1 = (conductance_MF(dF_1, T)/lF_1)*(frac) + (conductance_LF(dF_1, P1, P_pump1)/lF_1)*(1-frac);
				# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
					# .format(conductance_MF(dF_2, T)/lF_2, conductance_LF(dF_2, P2, P_pump2)/lF_2, C_2));
				# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
				# C_1 = conductance_LF(dF_1, P1, P_pump1)/lF_1;   # viscous flow
				# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
					# .format(x, P2, C_2))				
				regime_flag = 1;
	        #endif

			C_1 = conductance(dF_1, P1, P_pump1, L = lF_1);

	        #2. calculate Seff
			S_eff1 = 1/((1/S_1) + (1/C_1));

			#3. calculate chamber pressure
			P1 = (2*Q_0)/S_eff1;

	        #4. calculate pressure at pump inlet
			P_pump1 = P1 - (2*Q_0)/C_1;

	        #5. calculate error and update pumping speed
			err = abs((S_1 - np.interp(P_pump1, GV80[0], GV80[1], left = 0, right = 0))*100/S_1);
			S_1 = (S_1 + np.interp(P_pump1, GV80[0], GV80[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable
		prev_state_1 = (S_1, P1, P_pump1);

		# inner jacket pumping
		regime_1.append(regime_flag)
		mfp_1.append(meanFreePath(P1, T));
		P_1.append(P1);
		F_1.append(C_1);
		# print("{:.2f} L/s inner jacket conductance".format(C_1))
		# print("{:.2f} L/s (pfeiffer equation)".format(conductance_LF(dF_1, P1, P_pump1, L = lF_1)))
		# print("{:.2f} mbar -> {:.2f} mbar".format(P1, P_pump1))

		##############################################################################################
		# inner jacket to outer jacket
		# print(i)
		# print(P_1)
		# print(P_1)
		# print("\ninner jacket to outer...")

		delta_1 = rarefractionParameter(R_1/1E3, P_1[i]*100, v_0);
		W_1 = getWFactor(delta_1, LR_ratio_1);
		M_1 = massflow_lim(R_1/1E3, P_1[i]*100, v_0)*W_1;
		Q_1 = massFlow_to_throughPut(M_1,T);
		
		err = 0;
		first = True;
		regime_flag = 0;
		d = 4*2.54;		# viscous transition around mean free path ~= 4 inches
		low_d 		= d/3;
		high_d 		= 3*d;
		while(first or err > 2):
			if(first):
				S_2 = prev_state_2[0];
				P2 = prev_state_2[1];
				# if (P2 < 1E-5):
				# 	P2 = 1E-5
				#endif
				P_pump2 = prev_state_2[2];
				first = False;
	        #endif

	        # #1. calculate conductance
			# if(meanFreePath(P2, T) > 8*2.54):	
			# 	C_2 = conductance_MF(dF_2, T)/lF_2;	            # molecular flow
			# 	print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(P2, T), P2, C_2));
			# 	regime_flag = 0;
			# else:
			# 	C_2 = conductance_LF(dF_2, P2, P_pump2)/lF_2;   # viscous flow
			# 	print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(P2, T), P2, C_2))				
			# 	regime_flag = 1;
	        # # endif

	        #1. calculate conductance
			x = meanFreePath(P2, T);
			if(x > high_d):	
			# 	C_2 = conductance_MF(dF_2, T)/lF_2;	            # molecular flow
			# 	# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		# .format(x, P2, C_2));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
			# 	# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
			# 	# print("x... {:.2f}".format(x));
			# 	frac = (x - low_d)/(high_d - low_d);
			# 	C_2 = (conductance_MF(dF_2, T)/lF_2)*(frac) + (conductance_LF(dF_2, P2, P_pump2)/lF_2)*(1-frac);
			# 	# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
			# 		# .format(conductance_MF(dF_2, T)/lF_2, conductance_LF(dF_2, P2, P_pump2)/lF_2, C_2));
			# 	# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
			# 	C_2 = conductance_LF(dF_2, P2, P_pump2)/lF_2;   # viscous flow
			# 	# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		# .format(x, P2, C_2))				
				regime_flag = 1;
	        # endif

			C_2 = conductance(dF_2, P2, P_pump2, L = lF_2);

	        #2. calculate Seff
			S_eff2 = 1/((1/S_2) + (1/C_2));

			#3. calculate chamber pressure
			P2 = (2*Q_1)/S_eff2;

	        #4. calculate pressure at pump inlet
			P_pump2 = P2 - (2*Q_1)/C_2;

	        #5. calculate error and update pumping speed
			err = abs((S_2 - np.interp(P_pump2, EPX500[0], EPX500[1], left = 0, right = 0))*100/S_2);
			S_2 = (S_2 + np.interp(P_pump2, EPX500[0], EPX500[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable...
		prev_state_2 = (S_2, P2, P_pump2);

		# outer jacket pumping
		regime_2.append(regime_flag);
		mfp_2.append(meanFreePath(P2, T));
		P_2.append(P2);
		F_2.append(C_2);
		# print("{:.2f} L/s outer jacket conductance".format(C_2))


		##############################################################################################
		# outer jacket to Lesker Chamber
		# print("\nouter jacket to lesker...")
		delta_2 = rarefractionParameter(R_2/1E3, P_2[i]*100, v_0);
		W_2 = getWFactor(delta_2, LR_ratio_2);
		M_2 = massflow_lim(R_2/1E3, P_2[i]*100, v_0)*W_2;
		Q_2 = massFlow_to_throughPut(M_2,T);

		err = 0;
		first = True;
		regime_flag = 0;
		d = 16*2.54;
		low_d 		= d/3;
		high_d		= d*3;
		while(first or err > 2):
			if(first):
				S_3 = prev_state_3[0];
				P3 = prev_state_3[1];
				# if (P3 < 1E-5):
				# 	P3 = 1E-5
				#endif
				P_pump3 = prev_state_3[2];
				first = False;
	        #endif

	        #1. calculate conductance
			# if(meanFreePath(P3, T) > 45*2.54):		# viscous edge at mean free path > 45 inches
			# 	# print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P3, T), P3));
			# 	C_3 = conductance_MF(dF_3, T)/lF_3;	            # molecular flow
			# 	regime_flag = 0;
			# else:
			# 	# print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(P3, T), P3));
			# 	C_3 = conductance_LF(dF_3, P3, P_pump3)/lF_3;   # viscous flow
			# 	regime_flag = 1;	        
			# #endif

			x = meanFreePath(P3, T);
			if(x > high_d):	
			# 	C_3 = conductance_MF(dF_3, T)/lF_3;	            # molecular flow
			# 	# print("molecular flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		# .format(x, P3, C_3));
				regime_flag = 0;
			elif (x >= low_d and x <= high_d):
			# 	# print("low_d... {:.2f}, d... {:.2f}, high_d...{:.2f}".format(low_d, d, high_d))
			# 	# print("x... {:.2f}".format(x));
			# 	frac = (x - low_d)/(high_d - low_d);
			# 	C_3 = (conductance_MF(dF_3, T)/lF_3)*(frac) + (conductance_LF(dF_3, P3, P_pump3)/lF_3)*(1-frac);
			# 	# print("transistion region... molecular C: {:.2e}, viscous: {:.2e}, averaged.. {:.2e}"
			# 		# .format(conductance_MF(dF_3, T)/lF_3, conductance_LF(dF_3, P3, P_pump3)/lF_3, C_3));
			# 	# print("fraction MF... {:.3f}, fraction V...{:.3f}".format(frac, 1-frac));
				regime_flag	= 0;
			else:
			# 	C_3 = conductance_LF(dF_3, P3, P_pump3)/lF_3;   # viscous flow
			# 	# print("viscous flow... (path length... {:.2f} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		# .format(x, P3, C_3))				
				regime_flag = 1;
	        # endif
			C_3 = conductance(dF_3, P3, P_pump3, L = lF_3);


	        #2. calculate Seff
			S_eff3 = 1/((1/S_3) + (1/C_3));

			#3. calculate chamber pressure
			P3 = (Q_2)/S_eff3;

	        #4. calculate pressure at pump inlet
			P_pump3 = P3 - (Q_2)/C_3;

	        #5. calculate error and update pumping speed
			err = abs((S_3 - np.interp(P_pump3, EPX500[0], EPX500[1], left = 0, right = 0))*100/S_3);
			S_3 = (S_3 + np.interp(P_pump3, EPX500[0], EPX500[1], left = 0, right = 0))/2;
		#endwhile

		# update previous state variable
		prev_state_3	 = (S_3, P3, P_pump3);

		# Lesker Chamber pumping
		regime_3.append(regime_flag);
		mfp_3.append(meanFreePath(P3, T));
		P_3.append(P3);
		F_3.append(C_3);
		# print("{:.2f} L/s lesker chamber conductance".format(C_3))


		i +=1;
		# print("\n")
	#endfor

	# for i in range(len(P_0)):
	# 	print("({:.1f}, {:.3e}, {:d}, {:.3e})".format(P_0[i], P_1[i], regime_1[i], mfp_1[i]));
	#endfor
	# print(P_0)
	# print(P_1);
	# print(regime_1);
	# print(mfp_1);

	return((P_0, P_1, P_2, P_3), (regime_1, regime_2, regime_3), (mfp_1, mfp_2, mfp_3), (F_1, F_2, F_3));
#endfunction

def plotCellSim(P, R, M, C):
	'''
	Given the helium cell simulation outputs, plots the result.

		Parameters:
			P = (P_0, P_1, P_2, P_3)
				P_0 	(array(float)): helium cell pressures [mbar]
				P_1 	(array(float)): inner jacket pressures [mbar]
				P_2 	(array(float)): outer jacket pressures [mbar]
				P_3 	(array(float)): lesker chamber pressures [mbar]

			R = (regime_1, regime_2, regime_3)
				regime_1 (array(int)): regime trace for inner jacket(1 for viscous flow, 0 for molecular)
				regime_2 (array(int)): regime trace for outer jacket(1 for viscous flow, 0 for molecular)
				regime_3 (array(int)): regime trace for lesker chamber(1 for viscous flow, 0 for molecular)

			M = (mfp_1, mfp_2, mfp_3)
				mfp_1 (array(float)): mean free path for each helium cell pressure point for inner jacket [cm]
				mfp_2 (array(float)): mean free path for each helium cell pressure point for outer jacket [cm]
				mfp_3 (array(float)): mean free path for each helium cell pressure point for lesker chamber [cm]

			C = (F_1, F_2, F_3)
				F_1 (array(float)): effective conductance of inner jacket pumping system [L/s]
				F_2 (array(float)): effective conductance of outer jacket pumping system [L/s]
				F_3 (array(float)): effective conductance of lesker chamber pumping system [L/s]
		Returns:
			(fig, ax)
			fig 	(figure): figure reference for the summary helium cell sim plot
			ax 		(Axes): Axes reference for the summary helium cell sim plot

	'''
	fig, ax = plt.subplots(figsize = (10,4))

	# Pressure axis
	ax.set_yscale("log");
	# ax.set_xscale("log")
	ax.set_ylabel("Pressure [mbar]");
	ax.set_xlabel("Helium cell Pressure [mbar]");
	ax.set_title("3 Stage Helium Cell Pressures");
	plt.grid(True)

	# meaN free path axis
	ax2 = ax.twinx()
	# ax2.set_yscale("log");
	# ax2.set_ylabel("Mean Free Path [cm]");
	ax2.set_ylabel("Conductance [L/s]");


	# l1, = ax.plot(P[0], P[1], label = "Inner Jacket", color = "red");
	# l2, = ax.plot(P[0], P[2], label = "Outer Jacket", color = "blue");
	# l3, = ax.plot(P[0], P[3], label = "Lesker Chamber", color = "green");

	P1_v = ([],[]);
	P2_v = ([],[]);
	P3_v = ([],[]);

	
	P1_mf = ([],[]);
	P2_mf = ([],[]);
	P3_mf = ([],[]);


	for i in range(len(P[0])):
		#P1...
		if (R[0][i] == 0):
			P1_mf[0].append(P[0][i])
			P1_mf[1].append(P[1][i])
		else:
			P1_v[0].append(P[0][i])
			P1_v[1].append(P[1][i])
		#endif

		#P2...
		if (R[1][i] == 0):
			P2_mf[0].append(P[0][i])
			P2_mf[1].append(P[2][i])
		else:
			P2_v[0].append(P[0][i])
			P2_v[1].append(P[2][i])
		#endif

		#P3
		if (R[2][i] == 0):
			P3_mf[0].append(P[0][i])
			P3_mf[1].append(P[3][i])
		else:
			P3_v[0].append(P[0][i])
			P3_v[1].append(P[3][i])
		#endif
	#endfor

	l1_mf, 	= ax.plot(P1_mf[0], P1_mf[1], label = "Inner Jacket", color = "red");
	l1_v, 	= ax.plot(P1_v[0], P1_v[1], color = "red", linestyle = "--");

	l2_mf, 	= ax.plot(P2_mf[0], P2_mf[1], label = "Outer Jacket", color = "blue");
	l2_v, 	= ax.plot(P2_v[0], P2_v[1], color = "blue", linestyle = "--");

	l3_mf, 	= ax.plot(P3_mf[0], P3_mf[1], label = "Lesker Chamber", color = "green");
	l3_v, 	= ax.plot(P3_v[0], P3_v[1], color = "green", linestyle = "--");

	# mean free path lines...
	l4, = ax2.plot(P[0], C[0], color = "red", linestyle = "-.", alpha = .4)
	l5, = ax2.plot(P[0], C[1], color = "blue", linestyle = "-.", alpha = .4)
	l6, = ax2.plot(P[0], C[2], color = "green", linestyle = "-.", alpha = .4)

	ax.legend();

	# ax.legend(handles=[l1, l2,l3])
	# print(P0)
	# print(P1)
	# print(P2)
	# print(P3)
	# plt.grid(True)

	# plt.show(block = True);
	return (fig,ax);
#endfunction
