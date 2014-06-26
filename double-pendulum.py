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
import dp_draw

from math  import *
from numpy import *

from dp_lagrangian  import *
from dp_hamiltonian import *


##
# @brief prints usage instructions to stderr and exits
#
def print_usage():
	output  = "Usage: %s [OPTIONS]\n\n" % os.path.basename(__file__)
	output += "    -h, --help                    prints these instructions\n"
	output += "    -v, --verbose                 activates verbose mode\n"
	output += "    -H, --hamiltonian             runs simulation using Hamilton's equations\n"
	output += "    -g, --gravity=ACCEL           sets the gravitational acceleration\n"
	output += "    -s, --time-step=STEP          sets the simulation time step\n"
	output += "    -m, --mass=MASS1,MASS2        sets the mass of each bob\n"
	output += "    -t, --theta=THETA1,THETA2     sets the initial angle of each bob\n"
	output += "    -w, --omega=OMEGA1,OMEGA2     sets the initial angular velocity of each bob\n"
	output += "    -L, --rodlen=LEN1,LEN2        sets the rod length for each bob\n"
	output += "        --geometry=WIDTH,HEIGHT   sets the window dimensions\n\n"
	output += "Keyboard shortcuts:\n"
	output += "    Up Arrow      increases the time step\n"
	output += "    Down Arrow    decreases the time step\n"
	output += "    v             toggles verbose mode\n"
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
		opts, args = getopt.getopt(sys.argv[1:], "hvHm:t:w:L:s:g:",
			["help", "verbose", "mass=", "theta=", "omega=", "rodlen=",
			"time-step=", "gravity=", "geometry=", "hamiltonian"])
	except getopt.GetoptError as err:
		print_error(str(err))

	verbose = False

	# by default, use the Euler-Lagrange equations to simulate the system
	lagrangian = True

	for opt, arg in opts:
		if opt == "-v":
			verbose = True
		elif opt in ("-H", "--hamiltonian"):
			lagrangian = False
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
	if lagrangian:
		S = dp_lagrangian(g, m1, m2, t1, t2, w1, w2, L1, L2)
	else:
		S = dp_hamiltonian(g, m1, m2, t1, t2, w1, w2, L1, L2)

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
	while dp_draw.draw(S, window, Nx, Ny, dt):

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

		# check window events: quit, resize, key presses
		for event in pygame.event.get():
			if event.type == pygame.QUIT:
				return False
			elif event.type == pygame.VIDEORESIZE:
				(Nx,Ny) = event.size
				window = pygame.display.set_mode((Nx, Ny), pygame.RESIZABLE)
			elif event.type == pygame.KEYDOWN:
				if event.unicode == u'v':
					verbose = not verbose

		# the up and down arrow keys decrease and increase dt respectively
		pressed_keys = pygame.key.get_pressed()
		if pressed_keys[273]:
			dt *= 1.05
		if pressed_keys[274]:
			dt /= 1.05


if __name__ == '__main__':
	main()
