from dateutil import parser
import numpy as np
from timeutil import jd

SECtoDEG = 0.0002777777777777778

def mod_to_icrs(r, ut1):
    """
    Rotate coordinates from J2000 mean-equator of date to ICRS
    :param r: MOD position vector
    :param ut1: datetime of the vector
    :return: ICRS position vector
    """
    # Number of Julian centuries from the base epoch
    t = (jd(ut1) - 2451545.0) / 36525

    # These units are in arcseconds
    zeta = 2306.2181 * t + 0.30188 * t**2 + 0.017998 * t**3
    theta = 2004.3109 * t - 0.42665 * t**2 - 0.041833 * t**3
    # z = zeta + 0.79280 * t**2 + 0.000205 * t**3
    z = 2306.2181 * t + 1.09468 * t**2 + 0.018203 * t**3

    # Convert to radians
    zeta, theta, z = [np.deg2rad(x * SECtoDEG) for x in (zeta, theta, z)]

    # Build rotation matrix
    a0 = np.cos(theta) * np.cos(z) * np.cos(zeta) - np.sin(z) * np.sin(zeta)
    a1 = np.sin(z) * np.cos(theta) * np.cos(zeta) + np.sin(zeta) * np.cos(z)
    a2 = np.sin(theta) * np.cos(zeta)

    b0 = -np.sin(zeta) * np.cos(theta) * np.cos(z) - np.sin(z) * np.cos(zeta)
    b1 = -np.sin(z) * np.sin(zeta) * np.cos(theta) + np.cos(z) * np.cos(zeta)
    b2 = -np.sin(theta) * np.sin(zeta)

    c0 = -np.sin(theta) * np.cos(z)
    c1 = -np.sin(theta) * np.sin(z)
    c2 = np.cos(theta)

    rot = np.array([[a0, a1, a2],
                    [b0, b1, b2],
                    [c0, c1, c2]])

    return np.dot(rot, r)


if __name__ == "__main__":
    ut1 = parser.parse('April 2, 2006, 00:00 UTC')
    print(r(ut1))
