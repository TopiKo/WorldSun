from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# Hours for one roundrip around sun.
hourInSideralYear = 365.256*24
# Note earth rotation period w.r.t stars is not exactly 24h.
# 24h is the earth rotation period w.r.t sun
eartRotPeriod = 23. + 56./60 + 4./3600

init_printing(use_unicode=True)



# Define variables
lat, lon, earth_tilt = var('lat lon tilt', real = True) #lat_d/360*2*pi, lon_d/360*2*pi
rad_earth, dist_to_sun = var('R_earth dSun', real = True, positive = True)
x,y,z = symbols('x y z', real = True)
theta, phi = symbols('theta phi', real = True)
t = symbols('t', real = True, positive = True)



# Direction vector from sun to earth.
xDirEarth = cos(t/hourInSideralYear*2*pi)
yDirEarth = sin(t/hourInSideralYear*2*pi)
REarth = Matrix([[xDirEarth, yDirEarth, 0]])

# Latitude plane
theta_tilt = earth_tilt/360*2*pi
poleEarth = Matrix([[sin(theta_tilt), 0, cos(theta_tilt)]])
# Basis vectors in latitude plane
earthPlane1 = Matrix([[cos(theta_tilt), 0, -sin(theta_tilt)]])
earthPlane2 = Matrix([[0, 1, 0]])

# Expand observer in lat plane and pole vectors.
phiEarth = t/eartRotPeriod*2*pi
# Direction vector from earth center to observer
observerPlane = (cos(phiEarth)*earthPlane1 + sin(phiEarth)*earthPlane2)*cos(lat) \
                + poleEarth*sin(lat)
# vector from sun to observer
observer = REarth*dist_to_sun + observerPlane*rad_earth


sunDir = REarth

sunPlane1 = observerPlane.cross(sunDir)
sunPlane2 = observerPlane.cross(sunPlane1)

sunPlane1_perp = sunPlane1 - sunPlane1.dot(sunDir)*sunDir
sunPlane2_perp = sunPlane2 - sunPlane2.dot(sunDir)*sunDir

dA = sunPlane1.cross(sunPlane2).norm()
dA_perp = sunPlane1_perp.cross(sunPlane2_perp).norm()

frac = dA_perp/dA


def functions(lat_d = 0, lon_d = 0):
    # set of fixed parameters
    radiusE = 6300
    dSun = 10000


    params_set = [(earth_tilt, 23.),
                  (lat, lat_d/180*pi),
                  (lon, lon_d/180*pi),
                  (rad_earth, radiusE),
                  (dist_to_sun, .00000)]

    func_wrap = {}
    # Transform into numpy functions for speed..
    func_wrap['observer_f'] = utilities.lambdify(t, observer.subs(params_set), 'numpy')
    func_wrap['observerPlane_f'] = utilities.lambdify(t, observerPlane.subs(params_set)[:], "numpy")
    func_wrap['sunDir_f'] = utilities.lambdify(t, sunDir.subs(params_set)[:], "numpy")
    func_wrap['frac_f'] = utilities.lambdify(t, frac.subs(params_set), "numpy")
    func_wrap['poleEarth_f'] = utilities.lambdify(t, poleEarth.subs(params_set)[:], "numpy")
    func_wrap['earthPlane1_f'] = utilities.lambdify(t, earthPlane1.subs(params_set)[:], "numpy")
    func_wrap['earthPlane2_f'] = utilities.lambdify(t, earthPlane2.subs(params_set)[:], "numpy")
    func_wrap['REarth_f'] = utilities.lambdify(t, REarth.subs(params_set)[:], "numpy")

    return func_wrap

def sunShine(func_wrap, startDay, dayHours):

    plot_orbit = False
    day = 24.*startDay
    energy = 0
    sunProj_plot = []
    fracs, times = [], []


    for i, hour in enumerate(dayHours):
        t_s = hour + day
        sunProj = np.zeros((2,3))

        sky = np.array(func_wrap['observerPlane_f'](t_s))
        toNorthPole = (np.array(func_wrap['poleEarth_f'](t_s)) - sky)*np.linalg.norm(func_wrap['poleEarth_f'](t_s))
        north = (toNorthPole - toNorthPole.dot(sky)*sky)
        north /= np.linalg.norm(north)
        east = np.cross(north,sky)

        sunDirArr = np.array(func_wrap['sunDir_f'](t_s))
        sunX = -sunDirArr.dot(east)
        sunY = -sunDirArr.dot(north)
        sunZ = -sunDirArr.dot(sky)

        thetaSun = (np.pi/2 - np.arccos(sunZ))/np.pi*180

        if sunX != 0:
            angle = np.arctan(sunY/sunX) / (np.pi*2)*360

            if 0 < sunY:
                # the four cases of unit circle: note that north is 0 deg and east 90
                if 0 < sunX:
                    phiSun = 90. - angle
                elif sunX < 0:
                    phiSun = 270. - angle
                else: phiSun = 0
            else:
                if sunX < 0:
                    phiSun = 90 - angle + 180
                elif 0 < sunX:
                    phiSun = 90 - angle
                else: phiSun = 180
        else:
            if sunY < 0:
                phiSun = 180
            else:
                phiSun = 0

        if i != 0:
            size = 10
            if 0 < thetaSun:
                areaFrac = func_wrap['frac_f'](t_s)
                fracs.append(areaFrac)
                projPhiSun = 90-phiSun
                xSun = np.cos(thetaSun/180*np.pi)*np.cos(projPhiSun/180*np.pi)
                ySun = np.cos(thetaSun/180*np.pi)*np.sin(projPhiSun/180*np.pi)
                zSun = -np.sin(thetaSun/180*np.pi)
                point3 = np.array([xSun/(1-zSun), ySun/(1-zSun)])
                energy += areaFrac*(dayHours[i] - dayHours[i-1])
                sunProj_plot.append(point3)
                times.append(t_s)


    if len(sunProj_plot) != 0:
        return np.array(sunProj_plot), np.array(fracs), np.array(times), energy
    else:
        return None, None, None, 0
    #dayEnergies[k] = energy
