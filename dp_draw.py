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


import pygame

from math  import *
from numpy import *


##
# @brief draws the pendulum system on a window
# @param S the pendulum system class
# @param window the window where the pendulum system should be shown
# @param Nx window width (in pixels)
# @param Ny window height (in pixels)
# @param dt the simulation time step
# @return always True
#
def draw(S, window, Nx, Ny, dt):

	m1 = S.m1;  m2 = S.m2
	t1 = S.t1;  t2 = S.t2
	L1 = S.L1;  L2 = S.L2

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

	# write the time step value on the window
	myfont = pygame.font.SysFont("Arial", 15)
	label = myfont.render("dt = %.3g" % dt, 1, (128,128,128))
	window.blit(label, (10, 10))

	# update the screen
	pygame.display.flip()

	return True
