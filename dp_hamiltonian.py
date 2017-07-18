#!/usr/bin/env python3

"""Double pendulum simulator based on its Hamilton's equations."""

import math
import numpy


class DoublePendulumHamiltonian:

    def __init__(self, g, m1, m2, t1, t2, w1, w2, L1, L2):
        """
        Constructs a double pendulum simulator based on its
        Hamilton's equations. Bob #1 is the one attached to the fixed
        pivot.

        g - The gravitational acceleration.
        m1 - The mass of bob #1.
        m2 - The mass of bob #2.
        t1 - The initial angle of bob #1.
        t2 - The initial angle of bob #2.
        w1 - The initial angular velocity of bob #1.
        w2 - The initial angular velocity of bob #2.
        L1 - The length of the rod for bob #1.
        L2 - The length of the rod for bob #2.
        """

        self.g = g
        self.m1 = m1
        self.m2 = m2
        self.t1 = t1
        self.t2 = t2
        self.L1 = L1
        self.L2 = L2

        # compute the initial canonical momenta
        self.p1 = (m1 + m2) * (L1**2) * w1 + \
            m2 * L1 * L2 * w2 * math.cos(t1 - t2)
        self.p2 = m2 * (L2**2) * w2 + m2 * L1 * L2 * w1 * math.cos(t1 - t2)

    def potential_energy(self):
        """Computes the potential energy of the system."""

        m1 = self.m1
        t1 = self.t1
        L1 = self.L1
        m2 = self.m2
        t2 = self.t2
        L2 = self.L2

        g = self.g

        # compute the height of each bob
        y1 = -L1 * math.cos(t1)
        y2 = y1 - L2 * math.cos(t2)

        return m1 * g * y1 + m2 * g * y2

    def kinetic_energy(self):
        """Computes the kinetic energy of the system."""

        m1 = self.m1
        t1 = self.t1
        L1 = self.L1
        m2 = self.m2
        t2 = self.t2
        L2 = self.L2

        # compute the angular velocity of each bob
        (w1, w2) = self.omega()

        # compute the kinetic energy of each bob
        K1 = 0.5 * m1 * (L1 * w1)**2
        K2 = 0.5 * m2 * ((L1 * w1)**2 + (L2 * w2)**2 +
                         2 * L1 * L2 * w1 * w2 * math.cos(t1 - t2))

        return K1 + K2

    def mechanical_energy(self):
        """
        Computes the mechanical energy (total energy) of the
        system.
        """

        return self.kinetic_energy() + self.potential_energy()

    def omega(self):
        """
        Computes the angular velocities of the bobs and returns them
        as a tuple.
        """

        m1 = self.m1
        t1 = self.t1
        p1 = self.p1
        L1 = self.L1
        m2 = self.m2
        t2 = self.t2
        p2 = self.p2
        L2 = self.L2

        C0 = L1 * L2 * (m1 + m2 * math.sin(t1 - t2)**2)

        w1 = (L2 * p1 - L1 * p2 * math.cos(t1 - t2)) / (L1 * C0)
        w2 = (L1 * (m1 + m2) * p2 - L2 *
              m2 * p1 * math.cos(t1 - t2)) / (L2 * m2 * C0)

        return (w1, w2)

    def hamilton_rhs(self, t1, t2, p1, p2):
        """
        Computes the right-hand side of the Hamilton's equations for
        the double pendulum and returns it as an array.

        t1 - The angle of bob #1.
        t2 - The angle of bob #2.
        p1 - The canonical momentum of bob #1.
        p2 - The canonical momentum of bob #2.
        """

        m1 = self.m1
        L1 = self.L1
        m2 = self.m2
        L2 = self.L2

        g = self.g

        C0 = L1 * L2 * (m1 + m2 * math.sin(t1 - t2)**2)
        C1 = (p1 * p2 * math.sin(t1 - t2)) / C0
        C2 = (m2 * (L2 * p1)**2 + (m1 + m2) * (L1 * p2)**2 -
              2 * L1 * L2 * m2 * p1 * p2 * math.cos(t1 - t2)) * \
            math.sin(2 * (t1 - t2)) / (2 * C0**2)

        # F is the right-hand side of the Hamilton's equations
        F_t1 = (L2 * p1 - L1 * p2 * math.cos(t1 - t2)) / (L1 * C0)
        F_t2 = (L1 * (m1 + m2) * p2 - L2 *
                m2 * p1 * math.cos(t1 - t2)) / (L2 * m2 * C0)
        F_p1 = -(m1 + m2) * g * L1 * math.sin(t1) - C1 + C2
        F_p2 = -m2 * g * L2 * math.sin(t2) + C1 - C2

        return numpy.array([F_t1, F_t2, F_p1, F_p2])

    def time_step(self, dt):
        """
        Advances one time step using RK4 (classical Runge-Kutta
        method).
        """

        m1 = self.m1
        t1 = self.t1
        p1 = self.p1
        L1 = self.L1
        m2 = self.m2
        t2 = self.t2
        p2 = self.p2
        L2 = self.L2

        # y is an array with the canonical variables (angles + momenta)
        y = numpy.array([t1, t2, p1, p2])

        # compute the RK4 constants
        k1 = self.hamilton_rhs(*y)
        k2 = self.hamilton_rhs(*(y + dt * k1 / 2))
        k3 = self.hamilton_rhs(*(y + dt * k2 / 2))
        k4 = self.hamilton_rhs(*(y + dt * k3))

        # compute the RK4 right-hand side
        R = 1.0 / 6.0 * dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        # update the angles and momenta
        self.t1 += R[0]
        self.t2 += R[1]
        self.p1 += R[2]
        self.p2 += R[3]
