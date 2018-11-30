# Reference Ellipsoid (WSG84)

Reference ellipsoidal gravitational and rotational model, instantiated with 
WSG84 values for Earth.

Expected error of calculated values for gravity is less than 0.02% in the
surface normal component and 0.07% in the transverse component, where
percentages are given relative to the total magnitude of the computed gravity
vector (assuming correct implementation).

This should be evaluated relative to the fact the value of felt gravity differs
from its nominal value (9.80665 m/s^2) by extremes of +0.28% and -0.44% over
Earth's surface.

Basic coordinate transformations are included for convenience, allowing for full
account of non–inertial effects of the (geo)synchronous reference frame.

# Installation

## PyPi

`pip3 install --user rellipsoid`

## From source

Clone or download the source from [github](https://github.com/raabrp/rellipsoid)
Numpy and (for testing) pytest are required.

# Usage

* class `Planet` defines the model.
* `earth` instantiates it with WSG84 values.

## Methods

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

# Discussion of model

A planet's rotation deforms its shape and mass distribution (and thus 
experienced gravity close to its surface). Solving for the equilibrium between
gravity and centrifugal force (i.e. an equipotential surface) yields an
ellipsoid, which we use to approximate the surface of Earth. This approximation
is standard and is used nearly universally as a reference coordinate system 
(e.g. for GPS or geospacial data sets). The current standard reference ellipsoid
for Earth is part of the "WGS84" standard maintained by IERS.

## Unimplemented Corrections 

Empirical refinement of the model is possible with a spherical harmonic
expansions of gravitational potential (the gradient of which is local gravity). 

Such corrections are appropriate for highly accurate predictions of surface
gravity or orbital operations in the neighborhood of Earth. At extreme
distances however, the ellipsoidal model, or even a point–mass model is 
preferable.

Reference A. demonstrates calculations with the current IERS standard
spherical harmonic coefficients for Earth: the Earth Gravitational Model
1996 (EGM96).

More refined coefficients for Earth are provided by the NASA GRACE mission
(to be followed by GRACE-FO in 2019), which have allowed observation of temporal
variations.

ftp://ftp.csr.utexas.edu/pub/grace/GGM05

# Values used for Earth (WSG 84)

The current ellipsoidal model for Earth is the World Geodetic System 1984
(WSG 84): parameters which completely specify instantiation of the model.

Values for WSG 84 model obtained from Ref. A.:

```
6378137,       # a, semi-major axis (equatorial) [m]
298.257223563, # 1/f, inverse flatness [dimensionless]
7.292115e-5,   # omega, angular velocity [rad s^-1]
3.9860009e14   # GM, Universal grav. const. * mass of Earth [m^3 s^-2]
```

## Notes on Values Used

### omega

Treated as average contstant angular velocity for Earth, ignoring precession.
For very precise orbital calculations, such nuances may not be lightly ignored.

### GM

Different values for GM are used for different purposes

* old value 3.986005e14
* new value 3.986004418e14
* value without atmosphere 3.9860009e14  <-- We use this

In this package, we have used the value without atmosphere for near–surface
gravity estimations, justified by appeal to the shell theorem (approximately).

For the sake of completeness, orbital computations should use the updated
value, however GPS broadcasts assume that receivers continue to use the old
value and encode appropriately (no accuracy is lost while complicated
software updates to many deployed receivers are avoided).
 
## Evaluating Uncertainty

From the GRACE data, we note that gravitational potential can differ from
the predictions of the ellipsoidal model by 0.02% near the surface of
Earth, though such deviations are typically less than 0.005%. Ignoring the
transverse components of gravity caused by such irregularities, this
translates to errors of up to 0.02% in the magnitude of the calculated gravity
vector.

Accounting for 
[vertical deflection](https://en.wikipedia.org/wiki/Vertical_deflection), we 
note maximum deflections at Earth's surface of up to 100 arc-seconds (4.8e-4 
radians), yielding 0.00001% uncertainty in the vertical component (as a 
percentage of the total magnitude) and 0.048% in the transverse component. Note
that we use the geodetic definition of "vertical" (normal to the reference
ellipsoid).

Finally, in consideration of our choice of value for GM (which differs from the
updated value by subtracting the mass of the atmosphere), we note that the mass
of the atmosphere is ~0.0001% of that of the planet and consider that the force
of gravity linear in the planet's mass (this is approximately true: only for
pure gravitation, not the effects of centrifugal force).

We expect the ellipsoidal gravitational model, as implemented here for Earth,
to be accurate to +/- 0.02% in the vertical component and +/- 0.07% in the
transverse component (where percentages are given relative to the total
magnitude of the calculated gravity vector) everywhere on or near Earth's
surface.

All of this assumes of course that our neither our implementation nor our
analysis is flawed. This package it NOT intended for mission-critical
calculations or simulation. See license.txt for more information.
