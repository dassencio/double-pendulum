#!/usr/bin/python
# -*- coding: utf-8 -*-


################################################################################
#
#    Copyright (c) 2014, Diego Assencio (http://diego.assencio.com)
#    All rights reserved.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


from math  import *
from numpy import *


class dp_lagrangian:

	########################################################################
	#
	#	CLASS MEMBER VARIABLES
	#
	########################################################################

	## masses, angles, angular velocities and rod lengths
	m1 = 0.0;  m2 = 0.0
	t1 = 0.0;  t2 = 0.0
	w1 = 0.0;  w2 = 0.0
	L1 = 0.0;  L2 = 0.0

	## gravitational acceleration
	g = 0.0


	########################################################################
	#
	#	CLASS MEMBER FUNCTIONS
	#
	########################################################################

	##
	# @brief class constructor
	# @param g  gravitational acceleration
	# @param m1 mass of bob #1
	# @param m2 mass of bob #2
	# @param t1 initial angle of bob #1
	# @param t2 initial angle of bob #2
	# @param w1 initial angular velocity of bob #1
	# @param w2 initial angular velocity of bob #2
	# @note bob #1 is the one attached to the fixed pivot
	#
	def __init__(self, g, m1, m2, t1, t2, w1, w2, L1, L2):

		self.g  = g
		self.m1 = m1;  self.m2 = m2
		self.t1 = t1;  self.t2 = t2
		self.w1 = w1;  self.w2 = w2
		self.L1 = L1;  self.L2 = L2

	## computes the potential energy of the system
	def potential_energy(self):

		m1 = self.m1;  t1 = self.t1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  L2 = self.L2;

		g = self.g

		# compute the height of each bob
		y1 = -L1*cos(t1)
		y2 = y1 - L2*cos(t2)

		return m1*g*y1 + m2*g*y2

	## computes the kinetic energy of the system
	def kinetic_energy(self):

		m1 = self.m1;  t1 = self.t1;  w1 = self.w1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  w2 = self.w2;  L2 = self.L2;

		# compute the kinetic energy of each bob
		K1 = 0.5*m1*(L1*w1)**2
		K2 = 0.5*m2*((L1*w1)**2 + (L2*w2)**2 + 2*L1*L2*w1*w2*cos(t1-t2))

		return K1 + K2

	## computes the mechanical energy (total energy) of the system
	def mechanical_energy(self):
		return self.kinetic_energy() + self.potential_energy()

	##
	# @brief computes the right-hand side of the Euler-Lagrange equations
	#        for the pendulum system
	# @param t1 angle of bob #1
	# @param t2 angle of bob #2
	# @param w1 angular velocity of bob #1
	# @param w2 angular velocity of bob #2
	# @return the right-hand side of the Euler-Lagrange equations as an array
	#
	def lagrange_rhs(self, t1, t2, w1, w2):

		m1 = self.m1;  L1 = self.L1;
		m2 = self.m2;  L2 = self.L2;

		g = self.g

		a1 = (L2/L1)*(m2/(m1+m2))*cos(t1-t2)
		a2 = (L1/L2)*cos(t1-t2)

		f1 = -(L2/L1)*(m2/(m1+m2))*(w2**2)*sin(t1-t2) - (g/L1)*sin(t1)
		f2 = (L1/L2)*(w1**2)*sin(t1-t2) - (g/L2)*sin(t2)

		g1 = (f1 - a1*f2) / (1 - a1*a2)
		g2 = (f2 - a2*f1) / (1 - a1*a2)

		return array([w1, w2, g1, g2])

	## advances one time step using RK4 (classical Runge-Kutta method)
	def time_step(self, dt):

		m1 = self.m1;  t1 = self.t1;  w1 = self.w1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  w2 = self.w2;  L2 = self.L2;

		# y is an array with the generalized coordinates (angles +
		# angular velocities)
		y = array([t1, t2, w1, w2])

		# compute the RK4 constants
		k1 = self.lagrange_rhs(*y)
		k2 = self.lagrange_rhs(*(y + dt*k1/2))
		k3 = self.lagrange_rhs(*(y + dt*k2/2))
		k4 = self.lagrange_rhs(*(y + dt*k3))

		# compute the RK4 right-hand side
		R = 1.0/6.0 * dt * (k1 + 2.0*k2 + 2.0*k3 + k4)

		# update the angles and angular velocities
		self.t1 += R[0]
		self.t2 += R[1]
		self.w1 += R[2]
		self.w2 += R[3]
