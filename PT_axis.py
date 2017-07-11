''' 
This script computes Pressure (P), Tension (T) and Null (B) axes
using fault geometries: Strike, Dip and Rake angles from Moment tensor
inversions or First motion P-wave arrivals. The figure is plotted in the 
lower hemisphere using an equal-area projection. 

Esteban J. Chaves - Volcanological and Seismological Observatory of Costa Rica, Universidad Nacional (OVSICORI-UNA)-
Copyright(c) Esteban J. Chaves 2016-2017

This is code is distributed under the License GNU LGPL. 
Please cite using DOI: 10.5281/zenodo.825322

Contact: esteban.chaves.sibaja@una.cr
Website: https://github.com/echavess / http://www.ovsicori.una.ac.cr/index.php/ovsicori/personal
'''

import numpy as np

class Tensor_axes(object):

	def __init__(self, strike, dip, rake):

		'''Important: Strike, dip and rake must be Numpy N arrays event with only 1 element''' 

		self.strike = strike
		self.dip = dip
		self.rake = rake

	def  PTB(self):

		'''Lower hemisphere equal-area projection''' 
		projection = -1

		T = []
		P = []

		for strike, dip, rake in zip(self.strike, self.dip, self.rake):

			'''Computing Nomal (n) and Slip (u) vectors from Strike, Dip and Rake elements''' 
			n1 = -1*np.sin(dip*np.pi/180.0)*np.sin(strike*np.pi/180.0)
			n2 = np.sin(dip*np.pi/180.0)*np.cos(strike*np.pi/180.0)
			n3 = -1*np.cos(dip*np.pi/180.0)
			n = [n1, n2, n3]

			u1 = np.cos(rake*np.pi/180.0)*np.cos(strike*np.pi/180.0) + np.cos(dip*np.pi/180.0)*np.sin(rake*np.pi/180.0)*np.sin(strike*np.pi/180.0)
			u2 = np.cos(rake*np.pi/180.0)*np.sin(strike*np.pi/180.0) - np.cos(dip*np.pi/180.0)*np.sin(rake*np.pi/180.0)*np.cos(strike*np.pi/180.0)
			u3 = -1*np.sin(rake*np.pi/180.0)*np.sin(dip*np.pi/180.0)
			u = [u1, u2, u3]

			P_osa = np.subtract(n, u) / np.linalg.norm(np.subtract(n, u))
			T_osa = np.add(n, u) / np.linalg.norm(np.add(n, u))


			if P_osa[2] > 0:
				P_osa[0] = -1*P_osa[0]
				P_osa[1] = -1*P_osa[1]
				P_osa[2] = -1*P_osa[2]

			if T_osa[0] > 0:
				T_osa[0] = -1*T_osa[0]
				T_osa[1] = -1*T_osa[1]
				T_osa[2] = -1*T_osa[2]

			angle = (np.arctan(np.abs(P_osa[0]/P_osa[1]))*180.0)/np.pi

			if P_osa[0] > 0 and P_osa[1] > 0: 
				P_azimuth = angle

			elif P_osa[0] > 0 and P_osa[1] < 0:
				P_azimuth = 180 - angle

			elif P_osa[0] < 0 and P_osa[1] < 0:
				P_azimuth = angle + 180

			elif P_osa[0] < 0 and P_osa[1] > 0:
				P_azimuth = 360 - angle

			P_theta = (np.arccos(np.abs(P_osa[2]))*180.0)/np.pi


			angle_2 = (np.arctan(np.abs(T_osa[0]/T_osa[1]))*180.0)/np.pi

			if T_osa[0] > 0 and T_osa[1] > 0:
				T_azimuth = angle_2

			elif T_osa[0] > 0 and T_osa[1] < 0:
				T_azimuth = 180 - angle_2

			elif T_osa[0] < 0 and T_osa[1] < 0:
				T_azimuth = angle_2 + 180

			elif T_osa[0] < 0 and T_osa[1] > 0:
				T_azimuth = 360 - angle_2

			T_theta = (np.arccos(np.abs(T_osa[2]))*180.0)/np.pi

			P_x = np.sqrt(2.0)*projection*np.sin(P_theta*np.pi/360.0)*np.sin(P_azimuth*np.pi/180)
			P_y = np.sqrt(2.0)*projection*np.sin(P_theta*np.pi/360.0)*np.cos(P_azimuth*np.pi/180)

			P_axes = [P_x, P_y]
			P.append(P_axes)

			T_x = np.sqrt(2.0)*projection*np.sin(T_theta*np.pi/360)*np.sin(T_azimuth*np.pi/180)
			T_y = np.sqrt(2.0)*projection*np.sin(T_theta*np.pi/360)*np.cos(T_azimuth*np.pi/180)

			T_axes = [T_x, T_y]
			T.append(T_axes)

		print("Solving the Jacobian...")

		return P, T
