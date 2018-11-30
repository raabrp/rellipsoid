from pytest import approx
import numpy as np

from .rellipsoid import _rotmat, Planet, earth

def test_rotmat():
    '''
    unit test our axis-angle representation for rotations, since we use it
    throughout
    '''

    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])

    def r(a, b, c):
        '''rotate a by b (pi/2), compare to c'''
        assert _rotmat(b, np.pi/2).dot(a) == approx(c)

    # cyclic
    r(x, y, -z); r(y, z, -x); r(z, x, -y)

    # cyclic, negative angle
    r(x, -y, z); r(y, -z, x); r(z, -x, y)

    # cyclic, negative vector
    r(-x, y, z); r(-y, z, x);  r(-z, x, y)

    # cyclic, negative angle, negative vector
    r(-x, -y, -z); r(-y, -z, -x); r(-z, -x, -y)

    # anti-cyclic
    r(y, x, z); r(z, y, x); (x, z, y)

    # anti-cyclic, negative angle
    r(y, -x, -z); r(z, -y, -x); r(x, -z, -y)

    # anti-cyclic, negative vector
    r(-y, x, -z); r(-z, y, -x); r(-x, z, -y)

    # anti-cyclic, negative angle, negative vector
    r(-y, -x, z); r(-z, -y, x); r(-x, -z, y)


# nominal test values (mass accounts for atmosphere)
# no need for pytest fixtures for something so simple
n = Planet(6378137, 298.257223563, 7.292115e-5, 3986004.418e8)
n.gamma_p = 9.8321849378

def test_nominal_earth():
    '''compare against derived WSG84 values (Ref A.)'''

    def t(x, y):
        assert x == approx(y)

    t(n.b, 6356752.3142)
    t(n.e2, 6.69437999014e-3)
    t(n.e_p, 8.2094437949696e-2)
    t(n.E, 5.2185400842339e5)
    t(n.b / n.a, 0.996647189335)
    t(n.gamma_e, 9.7803253359)
    t(n.k, 0.00193185265241)
    t(n.m, 0.00344978650684)

def test_rad_curve_prime_vert():
    assert n._rad_curve_prime_vert(np.pi/2) == approx(6399593.6258)

def test_somigliana():
    def t(s2, g):
        assert n._somigliana(s2) == approx(g)

    t(0, n.gamma_e)
    t(1, n.gamma_p)

def test_free_air_gravity():
    def t(phi, h, g):
        assert n.get_free_air_gravity(phi, h) == approx(g)

    t(0, 0, n.gamma_e)
    t(np.pi/2, 0, n.gamma_p)
    t(-np.pi/2, 0, n.gamma_p)

def test_cartesian_geodetic():
    '''
    coordinate transformations
    Note the loss of accuracy in the inverse transform (geo -> cartesian)
    '''
    def t(phi, h, r, z):
        assert n._2d_cartesian_to_geodetic(r, z) == approx((phi, h))
        assert n._2d_geodetic_to_cartesian(phi, h) == approx((r, z), abs=1e-9)

    # equator
    t(0, 0, n.a, 0)
    t(0, 1000, n.a + 1000, 0)
    t(0, -1000, n.a - 1000, 0)

    # poles
    t(np.pi/2, 0, 0, n.b)
    t(-np.pi/2, 0, 0, -n.b)
    t(np.pi/2, 1000, 0, n.b + 1000)
    t(-np.pi/2, 1000, 0, -n.b - 1000)
    t(np.pi/2, -1000, 0, n.b - 1000)
    t(-np.pi/2, -1000, 0, -n.b + 1000)

def test_analytic_gravity():

    def t(phi, h, g):
        assert n.get_analytic_gravity(phi, h) == approx(g)

    t(0, 0, (-n.gamma_e, 0))
    t(np.pi/2, 0, (-n.gamma_p, 0))
    t(-np.pi/2, 0, (-n.gamma_p, 0))

def test_local_cartesian():
    '''
    More coordinate transforms. Again, geodetic -> cartesian is less accurate
    '''

    def t(phi, az, h):
        toloc, togeo = n.prep_local_cartesian(phi, az, h)
        def s(x, y, z, phi, lamda, h):
            # high absolute tolerance (2m)
            # for _2d_geodetic_to_cartesian
            assert toloc(phi, lamda, h) == approx((x, y, z), abs=1)

            p, l, j = togeo(x, y, z)
            assert p == approx(phi, abs=1e-6)

            # high absolute tolerance (2m)
            # for _2d_geodetic_to_cartesian
            assert j == approx(h, abs=2)

            # allow for degeneracy
            if not (p == approx(np.pi / 2) or p == approx(-np.pi / 2)):
                assert l == approx(lamda)
        return s

    s = t(0, 0, 0) # at equator looking north
    s(0, 0, 0, 0, 0, 0)                   # origin
    s(0, n.b, -n.a, np.pi/2, np.pi/2, 0)  # north pole
    s(0, -n.b, -n.a, -np.pi/2, 0, 0)      # south pole
    s(n.a, 0, -n.a, 0, np.pi/2, 0)        # east extreme
    s(-n.a, 0, -n.a, 0, -np.pi/2, 0)      # west extreme
    s(0, 0, -2*n.a, 0, np.pi, 0)          # opposite side

    # sanity check, if not numerical accuracy

    # at 40 latitude, y axis is west (x is north)
    p = 40 * np.pi / 180
    toloc, togeo = n.prep_local_cartesian(
        p, -np.pi/2, 0
    )
    phi, lamda, h = togeo(0, 0, 0) # at origin
    assert phi == approx(p)
    assert lamda == approx(0)
    assert h == approx(0)
    phi, lamda, h = togeo(0, 1e14, 0) # west (has radial component)
    assert phi < p
    assert lamda < 0
    assert h > 0
    phi, lamda, h = togeo(0, -1e14, 0) # east (has radial component)
    assert phi < p
    assert lamda > 0
    assert h > 0
    phi, lamda, h = togeo(0, 0, 1e14) # straight up
    assert h > 0
    assert phi == approx(p)
    assert lamda == approx(lamda)
    phi, lamda, h = togeo(0, 0, -2000) # straight down
    assert h < 0
    assert phi == approx(p)
    assert lamda == approx(lamda)
    phi, lamda, h = togeo(1e14, 0, 0) # north
    assert phi > p
    assert lamda == approx(lamda)
    assert h > 0
    phi, lamda, h = togeo(-1e14, 0, 0) # south
    assert phi < p
    assert lamda == approx(lamda)
    assert h > 0

    x, y, z = toloc(50 * np.pi / 180, 0, 0) # north
    assert x > 0
    assert z < 0
    assert y == approx(0, abs=1e-9) # again, geo -> cartesian
    x, y, z = toloc(30 * np.pi / 180, 0, 0) # south
    assert x < 0
    assert z < 0
    assert y == approx(0, abs=1e-9) # again, geo -> cartesian
    x, y, z = toloc(40 * np.pi / 180, -1, 0) # west
    assert x > 0
    assert z < 0
    assert y > 0
    x, y, z = toloc(40 * np.pi / 180, -4, 0) # east
    assert x > 0
    assert z < 0
    assert y < 0

def test_inertial():

    # at equator, y is north
    a, toN, toF = n.prep_local_cartesian_inertial(0, 0, 0)

    ax, ay, az = a(1 / n.omega) # medium interval of time
    assert ax < 0
    assert az < 0
    assert ay == approx(0, abs=1e-9)

    ax, ay, az = a(np.pi / n.omega) # half rotation of earth
    assert ax == approx(-n.a * np.pi)
    assert ay == approx(0, abs=1e-8)
    assert az == approx(-n.a * 2)

    nx, ny, nz = toN(np.pi / 2 / n.omega, [0, 0, 0]) # quarter rotation
    assert nz > 0
    assert nx == approx(-n.a)
    assert ny == approx(0, abs=1e-9)

    nx, ny, nz = toN(np.pi / n.omega, [0, 0, 0]) # half rotation
    assert nz == approx(-2 * n.a)
    assert nx == approx(-n.a * np.pi)
    assert ny == approx(0, abs=1e-8)

    # mid-atitude
    a, toN, toF = n.prep_local_cartesian_inertial(
        40 * np.pi / 180,
        0, 0
    )
    x, y, z = toN(10, (0, 1000, 0)) # flying north @ 100 m/s
    assert x > 0    # deflected right (decreasing radius)
    assert y < 1000 # surface is rotating around poll, viewed from lat > 0
    assert z > 0    # pulled upward

    x, y, z = toN(10, (1000, 0, 0)) # flying east @ 100 m/s
    assert x < 1000  # effect of rotation > slower velocity in inertial x component
    assert y < 0     # surface is rotating around poll, viewed from lat > 0
    assert z > 0     # pulled upward    
