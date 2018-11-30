#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
MIT License

Copyright (c) 2018 Reilly Raab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

###############################################################################

Reference ellipsoidal gravitational and rotational model, instantiated with
WSG84 values for Earth

Basic coordinate transformations are included for convenience, allowing for
full account of non–inertial effects of the (geo)synchronous reference frame.

# Usage:

* class `Planet` defines the model.
* `earth` instantiates it with WSG84 values.

## Methods:

    get_free_air_gravity: calculate magnitude of free air gravity (including
        centrifugal force) near the surface of the reference ellipsoid at
        specified latitude and height (geodetic)

    get_analytic_gravity: calculate gravity vector at a specified latitude
        and height (geodetic); contribution of centrifugal force is optional

    prep_local_cartesian: return transforms to and from a local Cartesian
        surface coordinate system and geodetic coordinates.

    prep_local_cartesian_inertial: return transformations between a local
        Cartesian surface coordinate system and an inertial frame
        coincident with it at time 0

# References

A. NIMA Technical Report TR8350.2, "Department of Defense World Geodetic System 
1984, Its Definition and Relationships With Local Geodetic Systems", Third 
Edition, Amendment 1, 3 January 2000: [link](http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf)

B. NIMA Technical Report TR8350.2, "Department of Defense World Geodetic System 
1984, Its Definition and Relationships With Local Geodetic Systems", Second
Edition: Chapters [3](http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%203.pdf) and
[4](http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%204.pdf)

C. Bernhard Hofmann-Wellenhof, Helmut Moritz; Physical Geodesy (2006)

D. Zhu, "Conversion of Earth-centered Earth-fixed coordinates to geodetic
   coordinates," Aerospace and Electronic Systems, IEEE Transactions on, vol. 30,
   pp. 957–961, 1994.

See readme.md for more information
'''

import numpy as np

def _rotmat(axis, theta):
    '''
    rotation matrix (around axis by theta)
    '''

    x, y, z = axis
    s, c = np.sin(theta), np.cos(theta)

    # rotation matrix from axis (n) and angle (theta)
    mat = np.array([
        [
            c + x ** 2 * (1 - c),
            x * y * (1 - c) - z * s,
            x * z * (1 - c) + y * s
        ],
        [
            y * x * (1 - c) + z * s,
            c + y ** 2 * (1 - c),
            y * z * (1 - c) - x * s
        ],
        [
            z * x * (1 - c) - y * s,
            z * y * (1 - c) + x * s,
            c + z ** 2 * (1 - c)
        ]
    ])

    return mat

###############################################################################

class Planet:
    '''
    Standard model and reference ellipsoidal coordinate system for the planet,
    used to approximate gravity and effects cause by a non-inertial (rotating)
    reference frame.

    Arguments:
        a:     semi-major axis (equatorial) [m]
        inv_f: inverse flatness [dimensionless]
        omega: angular velocity [rad s-1]
        gm:    Universal grav. const. * mass of planet [m^3 s^-2]


    Methods:

    get_free_air_gravity: calculate magnitude of free air gravity (including
        centrifugal force) near the surface of the reference ellipsoid at
        specified latitude and height (geodetic)

    get_analytic_gravity: calculate gravity vector at a specified latitude
        and height (geodetic); contribution of centrifugal force is optional

    prep_local_cartesian: return transforms to and from a local Cartesian
        surface coordinate system and geodetic coordinates.

    prep_local_cartesian_inertial: return transformations between a local
        Cartesian surface coordinate system and an inertial frame
        coincident with it at time 0
    '''
    def __init__(self, a, inv_f, omega, gm):

        f = 1 / inv_f

        self.a = a          # semi-major axis (equatorial radius) [m]
        self.f = f          # flattening []
        self.omega = omega  # angular velocity [rad s^-1]
        self.gm = gm        # Grav. const. * mass [m^3 s^-2]

        b = a * (1 - f)     # semi-minor axis [m]

        # first eccentricity squared
        e2 = f * (2 - f)
        #  = (a^2 - b^2) / a^2
        #  = E / a

        # second eccentricity
        e_p = np.sqrt(e2 / (1 - e2))
        #   = sqrt[ (a^2 - b^2) / b^2 ]
        #   = sqrt( e2 ) * (a / b)
        #   = E / b

        # linear eccentricity (half distance between foci)
        E = a * np.sqrt(e2)
        # = b * e_p
        # = sqrt[ a^2 - b^2 ]

        atan_ep = np.arctan(e_p) #  used twice below

        q_p = self._q_p(e_p)                                     # Ref. B(3-65)
        q0 = self._q(e_p)                                        # Ref. B(3-66)

        m = b * (omega * a) ** 2 / gm                            # Ref. B(4-10)

        # local gravity at pole [m s^-2]                         # Ref. B(3-63)
        gamma_e = gm / (a * b) * (1 - m - (m * e_p * q_p / (6 * q0)))

        # local gravity at equator [m s^-2]                      # Ref. B(3-64)
        gamma_p = gm / a ** 2 * (1 + m * e_p * q_p / (3 * q0))

        # constant for Somigliana equation
        k = (b * gamma_p) / (a * gamma_e) - 1

        # store derived values

        # purely geometric
        self.b = b
        self.e2 = e2
        self.e_p = e_p
        self.E = E
        self.q0 = q0
        self.q_p = q_p

        # depend on gm
        self.m = m
        self.gamma_e = gamma_e
        self.k = k

    @staticmethod
    def _q(x):
        '''
        A standard abbreviation in texts, where x has the form:
        (linear eccentricity / semi-minor axis) (dimensionless, in (0, inf))

        q = ((1 + 3 / x ** 2) * np.arctan(x) - (3 / x)) / 2

        The returned value is strictly positive, in (0, pi / 4)

        Ref. B(3-66); C(2-113)
        '''
        return ((1 + 3 / x ** 2) * np.arctan(x) - (3 / x)) / 2

    @staticmethod
    def _q_p(x):
        '''
        A standard abbreviation in texts, where x has the form:
        (linear eccentricity / semi-minor axis) (dimensionless, in (0, inf))

        q' = 3 * (1 + 1 / x ** 2) * (1 - np.arctan(x) / x) - 1

        The returned value is strictly positive, in (0, 2)

        Ref. B(3-65); C(2-133)
        '''
        return 3 * (1 + 1 / x ** 2) * (1 - np.arctan(x) / x) - 1

    def _rad_curve_prime_vert(self, phi):
        '''
        radius of curvature of the prime vertical [m], as a function of
        geodetic latitude [rad]
        '''

        a2 = self.a ** 2
        b2 = self.b ** 2

        s_phi, c_phi = np.sin(phi), np.cos(phi)

        return a2 / np.sqrt(a2 * c_phi ** 2 + b2 * s_phi ** 2)


    def _2d_geodetic_to_cartesian(self, phi, h):
        '''
        Convert geodetic to 2d Cartesian coordinates, ignoring longitude

        Note, our model is longitude-independent, thus we will only anticipate
        using r^2 = x^2 + y^2, but we retain this method for completeness.

        Arguments:

        phi: geodetic latitude [rad]
        h: geodetic height [m]

        Returns: cartesian coordinates (r, z)

        r: distance of point to rotational axis [m]
        z: distance of point to equatorial plane [m]
        '''
        n = self._rad_curve_prime_vert(phi)

        r = ((n + h) * np.cos(phi))
        z = (n * (self.b / self.a) ** 2 + h) * np.sin(phi)

        return r, z

    def _2d_cartesian_to_geodetic(self, r, z):
        '''
        Convert 2d Cartesian to geodetic coordinates, ignoring longitude

        Arguments:

        r: distance of point from rotational axis [m]
        z: distance of point from equatorial plane [m]

        Returns:

        phi: geodetic latitude [rad]
        h: geodetic height [m]
        '''
        # Uses method of Heikkinen
        # https://en.wikipedia.org/wiki/Geographic_coordinate_conversion

        r2 = r ** 2
        z2 = z ** 2
        e2 = self.e2
        e_p2 = self.e_p ** 2
        E2 = self.E ** 2
        F = 54 * self.b ** 2 * z2
        G = r2 + (1 - e2) * z2 - e2 * E2
        c = e2 ** 2 * F * r2 / G ** 3
        s = (1 + c + np.sqrt(c ** 2 + 2 * c)) ** (1 / 3)
        P = F / (3 * (s + 1 / s + 1) ** 2 * G ** 2)
        Q = np.sqrt(1 + 2 * e2 ** 2 * P)
        t1 = -(P * e2 * r) / (1 + Q)
        t2_1 = self.a ** 2 * (1 + 1 / Q) / 2
        t2_2 = P * (1 - e2) * z2 / (Q * (1 + Q))
        t2_3 = P * r2 / 2

        rad_val = t2_1 - t2_2 - t2_3
        if rad_val > 0:
            r0 = t1 + np.sqrt(rad_val)
        else:
            r0 = t1
        U = np.sqrt((r - e2 * r0) ** 2 + z2)
        V = np.sqrt((r - e2 * r0) ** 2 + (1 - e2) * z2)
        z0 = self.b ** 2 * z / (self.a * V)

        h = U * (1 - self.b ** 2 / (self.a * V))

        if r > 0:
            phi = np.arctan((z + e_p2 * z0) / r)
        else:
            phi = np.pi / 2 * np.sign(z)

        return phi, h

    def _2d_cartesian_to_harmonic(self, r, z):
        '''
        convert geodetic to ellipsoidal harmonic coordinates, ignoring
        longitude

        Arguments:

        r: distance of point from rotational axis [m]
        z: distance of point from equatorial plane [m]

        Returns:

        u: semi-minor axis of confocal ellipse passing through point [m]
        beta: reduced latitude [rad]
        '''

        r2 = r ** 2

        j = (r2 + z ** 2 - self.E ** 2)

        u = np.sqrt(j / 2 * (1 + np.sqrt(1 + 4 * (self.E * z / j) ** 2)))
        beta = np.arctan(z / u * np.sqrt((u ** 2 + self.E ** 2) / r2))

        return u, beta

    def _somigliana(self, s2):
        '''
        Somigliana equation

        Arguments:

        s2: sin^2(phi), where phi is geodetic latitude [rad]

        Returns:

        local apparent gravity [m s^-2] on the surface of the reference
        ellipsoid, accounting for centrifugal force.
        '''
        # Ref B(4-1)
        return self.gamma_e * (1 + self.k * s2) / np.sqrt(1 - (self.e2 * s2))

    def get_free_air_gravity(self, phi, h=0):
        '''
        Compute apparent gravity near the ellipsoidal surface, accounting for
        centrifugal force and height (free air correction), but not the
        Bouguer correction (assumes density of rock between point and
        reference ellipsoidal surface).

        Arguments:

        phi: geodetic latitude [rad]
        h: geodetic height [m]

        Returns:

        magnitude of local experienced gravity (and centrifugal force) [m s-2]
        '''

        s2 = np.sin(phi) ** 2
        gamma = self._somigliana(s2)

        # uses a Taylor expansion in height
        # eq. A(4-3) C(2-215)
        return gamma * (
            1 -
            2 * h * (1 + self.f + self.m - 2 * self.f * s2) / self.a +
            3 * (h / self.a) ** 2
        )

    def get_analytic_gravity(self, phi, h=0, centrifugal=True):
        '''
        Compute vector value of gravity, with optional centrifugal force.
        Applicable at all altitudes.

        Arguments:

        phi: geodetic latitude [rad]

        h: geodetic height [m]

        centrifugal: include centrifugal force (for geosynchronous reference
                     frames, such as the surface) [bool]

        Returns: v, n

        v: vertical component of (effective) gravitational field [m/s^2]
            (normal to ellipsoidal surface from which height is given)
        n: northward component of (effective) gravitational field [m/s^2]
            (relative to surface from which height is given)

        Intuition:

        As an extreme example, imagine a rapidly rotating, highly oblate
        planet. An observer looking "up" against the combined influence of
        gravity and centrifugal force at their location will observe an object
        apparently attracted to the poles by pure gravitation. The object, at a
        greater distance from the axis of rotation than the observer, also
        experiences greater centrifugal force than the surface-bound observer
        in any reference frame rotating with the planet, directed "up" (from
        the observer's perspective) and away from the nearest pole.

                                           Y
                                          ╱|
                                         ╱|
                                        ╱ |
                                    N  ╱ |
                              _ .. ```╱..X_
                            ╱   ╱   |╱  ╲   ╲
                           (   |   w->   |   )
                            ╲   ╲   |   ╱   ╱
                              ` ---...--- `
                                    S
        '''

        # coordinates
        r, z = self._2d_geodetic_to_cartesian(phi, h)   # cartesian
        u, beta = self._2d_cartesian_to_harmonic(r, z)  # harmonic
        s_beta = np.sin(beta) # dimensionless
        u2 = u ** 2           # [m2]

        # Ref A(4-10), C(2-128)
        # dimensionless > 0
        w = np.sqrt((u2 + (self.E * s_beta) ** 2) / (u2 + self.E ** 2))

        # Ref C(2-113)
        # dimensionless > 0
        q = self._q(self.E / u)

        # semimajor axis of confocal ellipse
        ac2 = (u ** 2 + self.E ** 2) # squared [m^2]
        ac = np.sqrt(ac2)           # [m]

        #
        # u-component of gravity (partial derivative of potential wrt
        # coordinate line normal to u-surface)
        #
        # negative value is "downward", as expected.
        # centrifugal component is "upward", as expected.
        #
        # Ref. C(2-132), A(4-5)
        #
        t1 = self.gm / ac2                                # [m s-2]
        t2_p1 = (self.omega * self.a) ** 2 * self.E / ac2 # [m s-2]
        t2_p2 = self.q_p / self.q0                        # dimensionless
        t2_p3 = np.sin(beta) ** 2 / 2 - 1/6               # dimensionless
        t2 = t2_p1 * t2_p2 * t2_p3
        t3 = self.omega ** 2 * u * np.cos(beta) ** 2      # [m s-2] > 0
        gamma_u = -(t1 + t2 - t3 * centrifugal) / w       # [m s-2]

        #
        # beta-component of gravity (as above, wrt beta.) The Jacobian for a
        # locally axis-aligned coordinate system is diagonal, so the only
        # transformation that takes place after partial differentiation of the
        # potential field is scaling.
        # Ref. C(2-131)
        #
        # Sign of pure gravitation is the same sign as beta (attraction towards
        # poles). Centrifugal force is opposite (away from poles).
        #
        # Ref. C(2-132) A(4-6)
        #
        t1_p1 = (self.omega * self.a) ** 2 / ac           # [m s-2] > 0
        t1_p2 = q / self.q0                               # dimensionless > 0
        t2 = self.omega ** 2 * ac                         # [m s-2] > 0
        p3 = np.sin(beta) * np.cos(beta)                  # dimensionless
        gamma_beta = -(-t1_p1 * t1_p2 + t2 * centrifugal) * p3 / w # [m s-2]

        #
        # geodetic latitude on the confocal ellipse containing the object.
        # gamma_u is normal to this confocal ellipse, the surface of which
        # forms angle phi_confocal with the equatorial plane.
        #
        # our expression for phi_confocal is positive. It is also strictly less
        # than phi in absolute value for h > 0 (implying a polar component of
        # gamma_u, as expected by the observation that confocal ellipses become
        # less eccentric - more circular - with greater size).
        #
        # Expression derived from Cartesian transform, setting h=0
        #
        phi_confocal = np.arctan(self.a * np.sqrt(ac2 - r ** 2) / (r * u))

        #
        # rotate gamma_u and gamma_beta into the surface geodetic coordinate
        # system
        #
        epsilon = phi - np.sign(phi) * phi_confocal
        gamma_h = gamma_u * np.cos(epsilon) - gamma_beta * np.sin(epsilon)
        gamma_phi = gamma_u * np.sin(epsilon) + gamma_beta * np.cos(epsilon)

        return gamma_h, gamma_phi

    def prep_local_cartesian(self, phi, az, h=0):
        '''
        Prepare transformations from geodetic coordinates to and from a
        surface-oriented cartesian coordinate system, where Z is normal to the
        surface, and Y points in the specified azimuthal direction

        Local Cartesian coordinates:

                               Ns  Y (On surface)
                                ╲  |
                            N-. p╲a|
                            ___`-.╲|______X (On surface)
                                  ╱|
                                 ╱ |
                                ╱  |
                          (Up) Z

        Ns: Cardinal North (on surface, XY plane)
        a:  azimuth (az) - angle between Y and Ns
        N:  Celestial North Pole (in the plane of Ns and Z)
        p:  latitude (phi) - angle between N and Ns

        Arguments:

        phi: geodetic latitude [rad]
        az: azimuth (clockwise from North on surface) [rad]
        h: geodetic height (default 0) [m]

        Returns: to_local, to_geodetic

        to_local: converts geodetic to local Cartesian coordinates
        to_geodetic: converts local Cartesian coordinates to geodetic

        All longitudes are relative to position at which the cartesian system
        was initialized
        '''

        sp, cp = np.sin(phi), np.cos(phi)
        sz, cz = np.sin(az), np.cos(az)

        # vector from geocenter to point on surface (assumed in XZ plane)
        r_com, z_com = self._2d_geodetic_to_cartesian(phi, h)

        def _local_to_geocentric(lx, ly, lz):
            '''
            converts local Cartesian to geocentric Cartesian, where origin of
            local coordinate system lies on reference meridian (0 longitude)
            '''
            # Rotate world around Z by (pi/2 - az) so X is South (surface), Y
            # is East (surface), and Z is normal to surface
            #
            #        Ns  Y (surf)                     Ns
            #         ╲  |                        N    |
            #     N-. p╲a|                         `. p|
            #     ___`-.╲|______X (surf)        _____`.|______Y (East)
            #           ╱|                            ╱|
            #          ╱ |                           ╱ |
            #         ╱  |                          ╱  |
            #   (Up) Z                             Z   X
            #
            # Rotate world around Y by (pi/2 - phi) so X parallel to Equatorial
            # plane, in the local meridian, Z is to Celestial North,
            # and Y is surface East
            #
            # Geocentric coordinate system, translated to point on surface at
            # reference meridian
            #
            #       X (Parallel to Equator)
            #        ╲         Z, N (North Celestial Pole)
            #         ╲  ╱ _.-`
            #     _____╲.-`p____Ns
            #          ╱
            #         ╱
            #        Y (East on Surface)
            #
            # p: (phi) latitude

            B = _rotmat([0, 0, 1], (np.pi / 2 - az))
            A = _rotmat([0, 1, 0], (np.pi / 2 - phi))

            rx, ry, rz = A.dot(B).dot([lx, ly, lz])

            # translate by vector from geocenter to local origin
            return np.array([rx + r_com, ry, rz + z_com])

        def _geocentric_to_local(x, y, z):
            # invert the above

            rx, ry, rz = x - r_com, y, z - z_com

            B = _rotmat([0, 1, 0], (phi - np.pi / 2))
            A = _rotmat([0, 0, 1], (az - np.pi / 2))

            return A.dot(B).dot([rx, ry, rz])

        def local_to_geodetic(x, y, z):
            '''
            convert local Cartesian to geodetic coordinates

            Arguments: x, y, z [m] (Cartesian vector in local coordinates)

            Returns: phi, lon, h

            phi: geodetic latitude [rad]
            lamda: relative longitude from local origin [rad]
            h: geodetic height [m]

            Local Cartesian coordinates:

                               Ns  Y (On surface)
                                ╲  |
                            N-. p╲a|
                            ___`-.╲|______X (On surface)
                                  ╱|
                                 ╱ |
                                ╱  |
                          (Up) Z

            Ns: Cardinal North (on surface, XY plane)
            a:  azimuth (az) - angle between Y and Ns
            N:  Celestial North Pole (in the plane of Ns and Z)
            p:  latitude (phi) - angle between N and Ns

            '''
            gx, gy, gz = _local_to_geocentric(x, y, z)

            r = np.sqrt(gx ** 2 + gy ** 2)

            phi, h = self._2d_cartesian_to_geodetic(r, gz)
            lamda = np.arctan2(gy, gx)

            return phi, lamda, h

        def geodetic_to_local(phi, lamda, h):
            '''
            convert geodetic to local Cartesian coordinates

            Arguments:

            phi: geodetic latitude [rad]
            lamda: relative longitude [rad]
            h: geodetic height [m]

            Returns: Cartesian vector in local coordinates [m, m, m]

            Local Cartesian coordinates:

                               Ns  Y (On surface)
                                ╲  |
                            N-. p╲a|
                            ___`-.╲|______X (On surface)
                                  ╱|
                                 ╱ |
                                ╱  |
                          (Up) Z

            Ns: Cardinal North (on surface, XY plane)
            a:  azimuth (az) - angle between Y and Ns
            N:  Celestial North Pole (in the plane of Ns and Z)
            p:  latitude (phi) - angle between N and Ns

            '''
            r, z = self._2d_geodetic_to_cartesian(phi, h)
            x = r * np.cos(lamda)
            y = r * np.sin(lamda)

            return _geocentric_to_local(x, y, z)

        return geodetic_to_local, local_to_geodetic

    def prep_local_cartesian_inertial(self, phi, az, h=0):
        '''
        For instances in which a local Cartesian coordinate system is
        appropriate, prepare transformations which will allow us to convert
        from an inertial reference to the noninertial one or vice-versa.

        We call the point on the surface of the rotating ellipsoid at which
        the origins of both references are coincident at time t=0 "A"

        The inertial frame is "F" and the noninertial frame is "N".

        Arguments:

        phi: geodetic latitude of A [rad]
        az: azimuth from of Y axis in local coordinate system [rad]
        h: geodetic height of A (default 0) [m]

        Returns: pos_A_in_F, transform_F_to_N, transform_N_to_F

        get_A_in_F: A function which returns the position of A as function of
                    time, in F

        transform_F_to_N: A function which transforms a vector in specified in
                          F at time t to a vector specified in N.

        transform_N_to_F: The inverse transformation of above.

        Explanation:

        Local Cartesian coordinates:

                   Ns  Y (On surface)
                    ╲  |
                N-. p╲a|                                ___ Direction of motion
                ___`-.╲|______X (On surface)          _.-`╱ relative to
                      ╱|                          _.-`      geocenter (vector
                     ╱ |                      _.-`          on surface of
                    ╱  |                   .-`              ellipsoid)
              (Up) Z

        Ns: Cardinal North (on surface, XY plane)
        a:  azimuth (az) - angle between Y and Ns
        N:  Celestial North Pole (in the plane of Ns and Z)
        p:  latitude (phi) - angle between N and Ns

        At time 0, we split this coordinate system into two coincident
        reference frames centered on A with zero relative velocity:

        * N: a non-inertial reference frame which rotates with the planet,
             maintaining A at its origin and aligning its Z axis to the surface
             normal and Y to the specified azimuthal direction for all time.
        * F: an inertial reference frame which moves at constant velocity,
             drifting off into space forever.

        '''

        # Solve for position of A in F

        sp, cp = np.sin(phi), np.cos(phi)
        sz, cz = np.sin(az), np.cos(az)

        n = self._rad_curve_prime_vert(phi)

        # distance from A to rotational axis of planet
        r = (n + h) * cp

        # Direction of North Celestial Pole in both F and N
        nx = -cp * sz
        ny = cp * cz
        nz = sp

        north = np.array([nx, ny, nz])

        def get_A_in_F(t):
            '''
            position of point A (origin of N), in frame F as function of time
            '''
            #
            # Solve in temporary coordinates with planet moving in +x direction
            # relative to F
            #
            #   |------ r w t ------|
            #
            #   F-------------- __.A0 -      Ao: point where frames were
            #                 A`    | |           coincident
            #                  ╲    | |       A: point on surface at later
            #                   f   | |          time (origin of N)
            #     Y              ╲wt| | r
            #     |               ╲ | |
            #     |___ X           ╲| |
            #    ╱             Axis + -
            #   Z
            #
            #  w: omega
            #  t: time
            #  r: distance from A to planet's axis of rotation
            #

            theta = self.omega * t

            x = r * (theta - np.sin(theta))
            y = r * (np.cos(theta) - 1)

            # Rotate to frame F
            #
            # Rotate world around X by (pi/2 - phi) so Z axis is normal to
            # surface at A0, Y is South, and X is West
            #
            #               Ns                 Ns  Y (On surface)
            #               |                   ╲  |
            #               |                    ╲a|
            # (West) X______|_______        ______╲|______X (On surface)
            #              ╱|                     ╱|
            #             ╱ |                    ╱ |
            #            ╱  |                   ╱  |
            #      (Up) Z   Y South       (Up) Z
            #
            # Rotate world around Z by (pi + az) so Y axis points to
            # given azimuth.
            #

            B = _rotmat([1, 0, 0], (np.pi / 2 - phi))
            A = _rotmat([0, 0, 1], (np.pi + az))

            return A.dot(B).dot([x, y, 0])

        def transform_F_to_N(t, v):
            '''
            Given a vector in inertial frame F and a time t,
            return the same vector in noninertial frame N

            Arguments:

            t: time [s]
            v: vector (x, y, z) [m] in frame F

            Returns:

            v as seen in N
            '''
            return _rotmat(north, -self.omega * t).dot(v - get_A_in_F(t))

        def transform_N_to_F(t, v):
            '''
            Given a vector in noninertial frame N and a time t,
            return the same vector in inertial frame F

            Arguments:

            t: time [s]
            v: vector (x, y, z) [m] in frame N

            Returns:

            v as seen in F
            '''
            return _rotmat(north, self.omega * t).dot(v) + get_A_in_F(t)

        return get_A_in_F, transform_F_to_N, transform_N_to_F

###############################################################################

earth = Planet(
    6378137,       # a, semi-major axis (equatorial) [m]
    298.257223563, # 1/f, inverse flatness []
    7.292115e-5,   # omega, angular velocity [rad s^-1]
    3.9860009e14   # GM, Universal grav. const. * mass of Earth [m^3 s^-2]
)
