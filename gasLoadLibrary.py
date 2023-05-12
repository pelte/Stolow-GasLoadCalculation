# Gas Load Library
# by: Patrick Elten
# 		started: 2023-02-13
# 		last updated:
# 
# This file will hold the functions facilitating the calculation of the gas load
# target pressures.
#
# initally created to validate the helium cell design of the HHG beamline.

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
ATOMIC_R	= 140					# Atomic Radius of He [pm]
COLL_A		= .21 					# Collision Cross Section [nm^2]
# VISCOSITY	= 1.825E-5				# Dynamic viscosity of air at 20 degrees C [Pa*s]


PI			= 3.1416
AVO_NUM 	= 6.0221367E23			# Avogadro's Number
K_CONST 	= 1.38064852E-23		# Boltzmann Constant [m^2kgs^-2K^-1]
R_CONST		= 83.1446261815324		# Universal Gas Constant [L*mbar*K^-1*mol^-1]


# given rarefraction parameter delta, duct length L, and aperture radius R
# return the dimensionless W factor. delta, L, R should be > 0.
# see Table 1 of "Rarefied gas flow through short tubes into vacuum" [DOI: 10.1116/1.2830639]
#
# This function DOES NOT EXTRAPOLATE IT ONLY INTERPOLATES. Given values outside the table range 
# will return a W factor at the nearest extreme.
def getWFactor(delta, in_r):
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
	ratioList = [0, 0.1, 0.5, 1, 5, 10];	# elements in this list are columns of Table 1
	deltaList = [0, .1, .5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000];	# elements in this list are rows of Table 1

	lines = [[],[],[],[],[],[]];
	print("\n")

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

# function returns the mass flow through an orifice at the free-molecular limit (delta = 0)
# M_0' [kg/s]
def massflow_lim(R, P_0, v_0):
	return m.sqrt(PI)*R*R*P_0/v_0;
#endfunction

# given the temperature returns the most probably molecular speed
# v_0 [m/s]
def molecularSpeed(T):
	return m.sqrt(2*K_CONST*T/ATOMIC_M_KG);
#endfunction;

# given aperture radius R [m], ref pressure P [Pa] and molecular speed [m/s] returns
# the refraction parameter \delta
def rarefractionParameter(R, P, v_0):
	return (R*P)/(MU_0*v_0);
#endfunction

def massFlow_to_throughPut(m, T):
	return m*R_CONST*T/MOLAR_M;
#endfunction

# returns the average speed of He at the given temprature
def averageSpeed(T):
	return m.sqrt(8*K_CONST*T/(PI*ATOMIC_M_KG));
#endfunction

# conductance [L/s] for one meter of a pipe of internal diameter D [m] experiencing molecular flow
def conductance_MF(D,T):
	return (2.6E-4*(averageSpeed(T)*100)*((D*100)**3))/100;
#endfunction

# conductance [L/s] for one meter of pipe of internal diameter D [m] experiencing laminar flow
# pressures given in mbar
def conductance_LF(D, P1, P2):
	return (1E3)*PI*(D**4)*(P1 + P2)*100/(256*MU_0*1);

# returns the required diameter of forelines given a target pumping speed and length
def dia_forelines(target, L, T):
	return (target / (averageSpeed(T)*100*2.6e-4/(100*L)))**(1/3)
#endfunction

# returns the beam radius [um] given a distance z [mm] from the focus
def beamRadius (z):
	f 		= 1 * 1E6;			# 1 m focal length [um]
	lmbd 	= 2;				# wavelength [um]
	half_div= .0012; 			# [rad]
	# M^2 = 1

	# waist radius [um]

	# using divergence:
	# w_0 = lmbd/(PI*half_div);

	# using beam diameter at lens:
	r_lens = (30/2) * 1E3;				# radius of beam on focusing lens [um]
	w_0 = (f*lmbd*1)/(PI*r_lens);
	# w_0 = (PI*d_lens*d_lens - m.sqrt(PI*PI*(d_lens**4) - 4*lmbd*lmbd*f*f*1))/(2*PI)
	print("w_0 {:.2f} um".format(w_0))

	# rayleigh range [um];
	z_R = (PI*(w_0**2))/lmbd;
	# print(w_0*)
	# print("z_R {:.2f} um".format(z_R))

	return w_0*m.sqrt(1+((z*1000/z_R)**2));
#endfunction

# returns the mean free path [cm] of He given T [K] and P [mbar]
def meanFreePath(T,P):
	return ((R_CONST*T/10)/(m.sqrt(2)*(COLL_A)*(1E-18)*AVO_NUM*(P*100)))*100;
#endfunction

# summary program
# performs the pressure calculation for 24 steps between the max and min cell pressures
def heliumCellSimulation(P_min, P_max, cell_length=5, num=24):
	T = 295;				# gas temperature [K] (assumed constant)
	tape_t = .0036*25.4;	# tape thicknmess [mm]

	P_0 = np.linspace(P_min, P_max, num=num, endpoint=True);
	P_1 = [];				# inner jacket pressure [mbar]
	P_2 = [];				# outer jacket pressure [mbar]
	P_3 = [];				# lesker chamber pressure [mbar]

	regime_1 = [];			# flow regime marker, 0 for molecular, 1 for viscous
	regime_2 = [];	
	regime_3 = [];

	mfp_1 = [];				# mean free path trace
	mfp_2 = [];
	mfp_3 = [];


	v_0 = molecularSpeed(T);

	# replace this with spot size calculator!
	R_0 = 2*beamRadius(cell_length/2)/ 1E3;
	print(R_0)
	L_0 = tape_t;
	LR_ratio_0 = L_0/R_0;
	S_initial_1 = 130;

	# this is a very optimistic conductance estimate for the purposes of this initial simulation
	# conductance in viscous or molecular flow regime?? add logic
	dF_1 = .1           # Diameter of pumping forelines [m]
	lF_1 = 10           # Length of pumping forelines [m]
	# C_1 = (conductance_MF(dF_1,T) / lF_1) # + iso 100 bellows inside chamber + gate valve
	# 									  # This should be viscous conductance!

	# S_eff1 = (1/(1/S_initial_1 + 1/C_1));

	R_1 = 2*beamRadius(.5*25.4)/ 1E3;
	print(R_1)
	L_1 = tape_t;
	LR_ratio_1 = L_1/R_1;
	S_initial_2 = 140;

	# this is a very optimistic conductance estimate for the purposes of this initial simulation
	# conductance in viscous or molecular flow regime?? add logic(?)
	dF_2 = .1           # Diameter of pumping forelines [m]
	lF_2 = 10           # Length of pumping forelines [m]
	# C_2 = (conductance_MF(dF_2,T) / lF_2) # + iso 100 bellows inside chamber + gate valve
	# S_eff2 = (1/(1/S_initial_2 + 1/C_2));

	R_2 = 3;				# conductance hole radius between outer jacket and chamber [mm]
	print(R_2);
	L_2 = 10*R_2;
	LR_ratio_2 = L_2/R_2;
	S_initial_3 = 140;

	# this is a very optimistic conductance estimate for the purposes of this initial simulation
	# conductance in viscous or molecular flow regime?? add logic
	dF_3 = .15          # Diameter of pumping forelines [m]
	lF_3 = 10           # Length of pumping forelines [m]
	# C_3 = (conductance_MF(dF_3,T) / lF_3) # + iso 100 bellows inside chamber + gate valve
	# S_eff3 = (1/(1/S_initial_3 + 1/C_3));

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
		print("\nhelium cell pressure... {:.2f} mbar".format(p));
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
			if(meanFreePath(T, P1) > 2*2.54):	# viscous edge at mean free path > 2 inches
				# print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(T,P1), P1));
				C_1 = conductance_MF(dF_1, T)/lF_1;	            # molecular flow
				regime_flag = 0;
			else:
				# print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(T,P1), P1));
				C_1 = conductance_LF(dF_1, P1, P_pump1)/lF_1;   # viscous flow
				regime_flag = 1;

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
		# P_1.append((2*Q_0)/S_eff1);
		regime_1.append(regime_flag)
		mfp_1.append(meanFreePath(T,P1));
		P_1.append(P1);

		##############################################################################################
		# inner jacket to outer jacket
		# print(i)
		# print(P_1)
		# print(P_1)
		print("\ninner jacket to outer...")

		delta_1 = rarefractionParameter(R_1/1E3, P_1[i]*100, v_0);
		W_1 = getWFactor(delta_1, LR_ratio_1);
		M_1 = massflow_lim(R_1/1E3, P_1[i]*100, v_0)*W_1;
		Q_1 = massFlow_to_throughPut(M_1,T);
		
		err = 0;
		first = True;
		regime_flag = 0;
		d = 2*4*2.54;		# viscous transition around mean free path ~= 4 inches
		low_d 		= d/5;
		high_d 		= 5*d;
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
			# if(meanFreePath(T, P2) > 8*2.54):	
			# 	C_2 = conductance_MF(dF_2, T)/lF_2;	            # molecular flow
			# 	print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(T,P2), P2, C_2));
			# 	regime_flag = 0;
			# else:
			# 	C_2 = conductance_LF(dF_2, P2, P_pump2)/lF_2;   # viscous flow
			# 	print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar, C = {:.2e})"
			# 		.format(meanFreePath(T,P2), P2, C_2))				
			# 	regime_flag = 1;
	        # # endif

	        #1. calculate conductance
			x = meanFreePath(T,P2);
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
		# P_2.append((2*Q_1)/S_eff2);
		print("final pressure... {:.2e}".format(P2))
		regime_2.append(regime_flag);
		mfp_2.append(meanFreePath(T,P2));
		P_2.append(P2);

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
		d = 2*16*2.54;
		low_d 		= d/10;
		high_d		= .001;
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
			# if(meanFreePath(T, P3) > 45*2.54):		# viscous edge at mean free path > 45 inches
			# 	# print("molecular flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(T,P3), P3));
			# 	C_3 = conductance_MF(dF_3, T)/lF_3;	            # molecular flow
			# 	regime_flag = 0;
			# else:
			# 	# print("viscous flow... (path length... {:.2e} cm, P = {:.2e} mbar)".format(meanFreePath(T,P3), P3));
			# 	C_3 = conductance_LF(dF_3, P3, P_pump3)/lF_3;   # viscous flow
			# 	regime_flag = 1;	        
			# #endif

			x = meanFreePath(T,P3);
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
		# P_3.append((Q_2)/S_eff3);
		regime_3.append(regime_flag);
		mfp_3.append(meanFreePath(T,P3));
		P_3.append(P3);

		i +=1;
	#endfor

	# for i in range(len(P_0)):
	# 	print("({:.1f}, {:.3e}, {:d}, {:.3e})".format(P_0[i], P_1[i], regime_1[i], mfp_1[i]));
	#endfor
	# print(P_0)
	# print(P_1);
	# print(regime_1);
	# print(mfp_1);

	return((P_0, P_1, P_2, P_3), (regime_1, regime_2, regime_3), (mfp_1, mfp_2, mfp_3));
#endfunction

# plot the results of the helium cell simulation
def plotCellSim(P, R, M):
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
	ax2.set_yscale("log");
	ax2.set_ylabel("Mean Free Path [cm]");

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
	l4, = ax2.plot(P[0], M[0], color = "red", linestyle = "-.", alpha = .3)
	l5, = ax2.plot(P[0], M[1], color = "blue", linestyle = "-.", alpha = .3)
	l6, = ax2.plot(P[0], M[2], color = "green", linestyle = "-.", alpha = .3)

	ax.legend();

	# ax.legend(handles=[l1, l2,l3])
	# print(P0)
	# print(P1)
	# print(P2)
	# print(P3)
	# plt.grid(True)

	plt.show(block = True);
	return (fig,ax);
#endfunction