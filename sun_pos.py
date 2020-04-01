from dateutil import parser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local utilities
from timeutil import jd
from coord_trans import mod_to_icrs

# Astropy imports
from astropy.time import Time
from astropy import coordinates

AUtoKM = 149597870.7    # 1 AU in kilometers


def pos(ut1, to_icrs=True):
    """
    Return the Sun position vector in J2000 Mean-equator of date frame
    Algorithm 29 from Vallado, 3rd edition (also Wikipedia, also Montenbruck)
    :param ut1: datetime
    :return: Sun position vector in an MOD (mean-equator of date) frame in km
    """
    # Julian centuries
    t = (jd(ut1) - 2451545.0) / 36525

    # Mean longitude of the Sun in degrees (in the range 0->360)
    lm = (280.460 + 36000.771 * t) % 360

    # Mean anomaly of the sun in degrees (in the range 0->360)
    m = (357.5277233 + 35999.05034 * t) % 360

    # anomaly of the ecliptic
    rad_m = np.deg2rad(m)
    le = lm + 1.914666471 * np.sin(rad_m) + 0.019994643 * np.sin(2 * rad_m)

    # Obliquitiy of the ecliptic
    ep = (23.439297 - 0.0130042 * t) % 360
    rad_ep = np.deg2rad(ep)

    # Distance to the sun (in AU)
    r_mag = 1.000140612 - 0.016708617 * np.cos(rad_m) - \
            0.000139589 * np.cos(2 * rad_m)

    # Vector to the sun (in AU)
    rad_le = np.deg2rad(le)
    r_vec = np.array([r_mag * np.cos(rad_le),
                      r_mag * np.cos(rad_ep) * np.sin(rad_le),
                      r_mag * np.sin(rad_ep) * np.sin(rad_le)]
                     )
    r_vec = r_vec * AUtoKM

    # Rotate to ICRS if required
    if to_icrs:
        return mod_to_icrs(r_vec, ut1)
    else:
        return r_vec


def astropy_validation(ut1):
    s = coordinates.get_sun(Time(ut1))
    return s.cartesian.xyz.to('km').value


def vector_angle(r1, r2):
    dot = np.dot(r1, r2)
    r1_mag, r2_mag = np.linalg.norm([r1, r2], axis=1)
    theta = np.arccos(dot / (r1_mag * r2_mag))

    return np.rad2deg(theta)


if __name__ == "__main__":
    # Print out values from the book
    ut1 = parser.parse('April 2, 2006, 00:00 UTC')
    r = pos(ut1)
    print(r)
    print(astropy_validation(ut1))

    # Validate against astropy and SPICE
    times = pd.date_range(start='1/1/2020', end='1/1/2030', freq='1D',
                          tz='UTC')
    df = pd.DataFrame(index=times)

    # Shove the different results into the dataframe
    r = list(zip(*[pos(x) for x in times]))
    df['vallado_x'] = r[0]
    df['vallado_y'] = r[1]
    df['vallado_z'] = r[2]

    r = list(zip(*[astropy_validation(x) for x in times]))
    df['astropy_x'] = r[0]
    df['astropy_y'] = r[1]
    df['astropy_z'] = r[2]

    # Grab SPICE data from:
    # https://wgc.jpl.nasa.gov:8443/webgeocalc/#StateVector
    # J2000, 1/1/2020-1/1/2021

    # Data starts on line 15, last 4 lines are parameters
    spice_df = pd.read_csv('WGC_StateVector_20200401122555.csv',
                           skiprows=15, skipfooter=4, parse_dates=True,
                           index_col='UTC calendar date', engine='python')

    # Join the data and get the error of the  Vallado and astropy vectors
    jdf = df.join(spice_df)
    jdf['vallado_error'] = jdf.apply(lambda x: vector_angle([x.vallado_x, x.vallado_y, x.vallado_z],
                                                            [x['X (km)'], x['Y (km)'], x['Z (km)']]),
                                     axis=1)

    jdf['astropy_error'] = jdf.apply(lambda x: vector_angle([x.astropy_x, x.astropy_y, x.astropy_z],
                                                            [x['X (km)'], x['Y (km)'], x['Z (km)']]),
                                     axis=1)

    # Just plot the errors
    fig, ax = plt.subplots()
    jdf[['vallado_error', 'astropy_error']].plot(ax=ax, alpha=0.5)
    ax.set_ylabel('Degrees off of SPICE ephemeris')
    fig.savefig('error.png')

