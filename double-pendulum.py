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


import os
import sys
import getopt
import pygame

from math  import *
from numpy import *


class pendulum_system:

	########################################################################
	#
	#	CLASS MEMBER VARIABLES
	#
	########################################################################

	## masses, angles, canonical momenta and rod lengths
	m1 = 0.0;  m2 = 0.0
	t1 = 0.0;  t2 = 0.0
	p1 = 0.0;  p2 = 0.0
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
		self.L1 = L1;  self.L2 = L2

		# compute the initial canonical momenta
		self.p1 = (m1+m2)*(L1**2)*w1 + m2*L1*L2*w2*cos(t1-t2)
		self.p2 = m2*(L2**2)*w2 + m2*L1*L2*w1*cos(t1-t2)

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

		m1 = self.m1;  t1 = self.t1;  p1 = self.p1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  p2 = self.p2;  L2 = self.L2;

		# compute the angular velocity of each bob
		(w1,w2) = self.omega()

		# compute the kinetic energy of each bob
		K1 = 0.5*m1*(L1*w1)**2
		K2 = 0.5*m2*((L1*w1)**2 + (L2*w2)**2 + 2*L1*L2*w1*w2*cos(t1-t2))

		return K1 + K2

	## computes the mechanical energy (total energy) of the system
	def mechanical_energy(self):
		return self.kinetic_energy() + self.potential_energy()

	## computes the angular velocity of each bob and returns them as a tuple
	def omega(self):

		m1 = self.m1;  t1 = self.t1;  p1 = self.p1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  p2 = self.p2;  L2 = self.L2;

		C0 = L1*L2*(m1 + m2*sin(t1-t2)**2)

		w1 = (L2*p1 - L1*p2*cos(t1-t2)) / (L1*C0)
		w2 = (L1*(m1+m2)*p2 - L2*m2*p1*cos(t1-t2)) / (L2*m2*C0)

		return (w1,w2)

	##
	# @brief computes the right-hand side of the Hamilton's equations for
	#        the pendulum system
	# @param t1 angle of bob #1
	# @param t2 angle of bob #2
	# @param p1 canonical momentum of bob #1
	# @param p2 canonical momentum of bob #2
	# @return the right-hand side of the Hamilton's equations as an array
	#
	def hamilton_rhs(self, t1, t2, p1, p2):

		m1 = self.m1;  L1 = self.L1;
		m2 = self.m2;  L2 = self.L2;

		g = self.g

		C0 = L1*L2*(m1 + m2*sin(t1-t2)**2)
		C1 = (p1*p2*sin(t1-t2)) / C0
		C2 = (m2*(L2*p1)**2 + (m1+m2)*(L1*p2)**2 - L1*L2*m2*p1*p2*cos(t1-t2)) * \
		     sin(2*(t1-t2)) / (2*C0**2)

		# F is the right-hand side of the Hamilton's equations
		F_t1 = (L2*p1 - L1*p2*cos(t1-t2)) / (L1*C0)
		F_t2 = (L1*(m1+m2)*p2 - L2*m2*p1*cos(t1-t2)) / (L2*m2*C0)
		F_p1 = -(m1 + m2)*g*L1*sin(t1) - C1 + C2
		F_p2 = -m2*g*L2*sin(t2) + C1 - C2

		return array([F_t1, F_t2, F_p1, F_p2])

	## advances one time step using RK4 (classical Runge-Kutta method)
	def time_step(self, dt):

		m1 = self.m1;  t1 = self.t1;  p1 = self.p1;  L1 = self.L1;
		m2 = self.m2;  t2 = self.t2;  p2 = self.p2;  L2 = self.L2;

		# y is an array with the canonical variables (angles + momenta)
		y = array([t1, t2, p1, p2])

		# compute the RK4 constants
		k1 = self.hamilton_rhs(*y)
		k2 = self.hamilton_rhs(*(y + dt*k1/2))
		k3 = self.hamilton_rhs(*(y + dt*k2/2))
		k4 = self.hamilton_rhs(*(y + dt*k3))

		# compute the RK4 right-hand side
		R = 1.0/6.0 * dt * (k1 + 2.0*k2 + 2.0*k3 + k4)

		# update the angles and momenta
		self.t1 += R[0]
		self.t2 += R[1]
		self.p1 += R[2]
		self.p2 += R[3]

	##
	# @brief draws the pendulum system on a window
	# @param window the window where the pendulum system should be shown
	# @param Nx window width (in pixels)
	# @param Ny window height (in pixels)
	# @return always True
	#
	def draw(self, window, Nx, Ny):

		m1 = self.m1;  m2 = self.m2
		t1 = self.t1;  t2 = self.t2
		L1 = self.L1;  L2 = self.L2

		# radius (in pixels) of each bob (min/max: 3/12 pixels)
		R1 = max(3, int( 12 * (m1 / (m1 + m2)) ))
		R2 = max(3, int( 12 * (m2 / (m1 + m2)) ))

		# length (in pixels) of each rod
		P1 = 0.85 * min(Nx/2,Ny/2) * (L1 / (L1 + L2))
		P2 = 0.85 * min(Nx/2,Ny/2) * (L2 / (L1 + L2))

		# positions (in (pixels,pixels)) of each bob
		X0 = array([Nx/2,Ny/2])
		X1 = X0 + array([int(P1*sin(t1)), int(P1*cos(t1))])
		X2 = X1 + array([int(P2*sin(t2)), int(P2*cos(t2))])

		# color: rods and bobs
		color_L1 = (255, 255, 255)
		color_L2 = (128, 128, 128)
		color_m1 = (255, 0  , 0  )
		color_m2 = (0  , 0  , 255)

		# clear the window
		window.fill((0,0,0))

		# draw the rods and the bobs
		pygame.draw.line(window, color_L1, X0, X1, 3)
		pygame.draw.line(window, color_L2, X1, X2, 3)
		pygame.draw.circle(window, color_m1, X1, R1)
		pygame.draw.circle(window, color_m2, X2, R2)

		# update the screen
		pygame.display.flip()

		return True


##
# @brief prints usage instructions to stderr and exits
#
def print_usage():
	output  = "Usage: %s [OPTIONS]\n\n" % os.path.basename(__file__)
	output += "    -h, --help                   prints these instructions\n"
	output += "    -v, --verbose                activates verbose mode\n"
	output += "    -g, --gravity=ACCEL          sets the gravitational acceleration\n"
	output += "    -s, --time-step=STEP         sets the simulation time step\n"
	output += "    -m, --mass=MASS1,MASS2       sets the mass of each bob\n"
	output += "    -t, --theta=THETA1,THETA2    sets the initial angle of each bob\n"
	output += "    -w, --omega=OMEGA1,OMEGA2    sets the initial angular velocity of each bob\n"
	output += "    -L, --rodlen=LEN1,LEN2       sets the rod length for each bob\n"
	output += "        --geometry=WIDTH,HEIGHT   sets the window dimensions\n\n"
	sys.stderr.write(output)
	sys.exit(0)


##
# @brief prints an error message and exits with an error code (1)
#
def print_error(errmsg):
	sys.stderr.write("Error: %s\n" % errmsg)
	sys.exit(1)


def main():

	# default simulation parameter values
	g  = 10
	dt = 0.01
	m1 = m2 = 1.0
	t1 = t2 = 0.5
	w1 = w2 = 0.0
	L1 = L2 = 1.0

	# default window size parameters
	Nx = Ny = 500

	#
	# process the input parameters
	#
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hvm:t:w:L:s:g:",
			["help", "verbose", "mass=", "theta=", "omega=", "rodlen=",
			"time-step=", "gravity=", "geometry="])
	except getopt.GetoptError as err:
		print_error(str(err))

	verbose = False

	for opt, arg in opts:
		if opt == "-v":
			verbose = True
		elif opt in ("-h", "--help"):
			print_usage()
		elif opt in ("-m", "--mass"):
			try:
				(m1,m2) = map(lambda x: float(x), tuple(arg.split(',')))
				if m1 <= 0 or m2 <= 0:
					raise ValueError
			except:
				print_error("invalid mass values (masses must be positive)")
		elif opt in ("-t", "--theta"):
			try:
				(t1,t2) = map(lambda x: float(x), tuple(arg.split(',')))
			except:
				print_error("invalid initial angle values")
		elif opt in ("-w", "--omega"):
			try:
				(w1,w2) = map(lambda x: float(x), tuple(arg.split(',')))
			except:
				print_error("invalid initial angular velocity values")
		elif opt in ("-L", "--rodlen"):
			try:
				(L1,L2) = map(lambda x: float(x), tuple(arg.split(',')))
				if L1 <= 0 or L2 <= 0:
					raise ValueError
			except:
				print_error("invalid rod length values (rod lengths must be positive)")
		elif opt in ("-s", "--time-step"):
			try:
				dt = float(arg)
				if dt <= 0:
					raise ValueError
			except:
				print_error("invalid time step value (time step must be positive)")
		elif opt in ("-g", "--gravity"):
			try:
				g = float(arg)
			except:
				print_error("invalid gravitational acceleration value")
		elif opt in ("--geometry"):
			try:
				(Nx,Ny) = map(lambda x: int(x), tuple(arg.split(',')))
				if Nx <= 0 or Ny <= 0:
					raise ValueError
			except:
				print_error("invalid window dimensions (dimensions must be positive)")
		else:
			print_usage()

	# initialize the pendulum system
	S = pendulum_system(g, m1, m2, t1, t2, w1, w2, L1, L2)

	# E0 = initial mechanical energy of the system
	E0 = S.mechanical_energy()

	step = 0

	# maximum energy change (compared to E0): if too large, the simulation is bad
	max_dE = 0

	pygame.init()

	clock = pygame.time.Clock()

	# create the output window
	window = pygame.display.set_mode((Nx, Ny), pygame.RESIZABLE)

	pygame.display.set_caption("double pendulum")

	# keep running the simulation until the user closes the window
	while S.draw(window, Nx, Ny):

		# limit the while loop to a max of 25 times per second.
		clock.tick(25)

		# advance one time step
		S.time_step(dt)

		if verbose:
			Et = S.mechanical_energy()
			max_dE = max(abs(Et - E0), max_dE)
			print "[%u] t = %f   Et = %f   E0 = %f   |Et - E0| = %f   |Et - E0|_max = %f" % (
				step, step*dt, Et, E0, abs(Et-E0), max_dE
			)

		step += 1

		# check window events: quit and resize
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				return False
			if event.type == pygame.VIDEORESIZE:
				(Nx,Ny) = event.size
				window = pygame.display.set_mode((Nx, Ny), pygame.RESIZABLE)


if __name__ == '__main__':
	main()
