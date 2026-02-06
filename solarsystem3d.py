import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
from matplotlib.animation import FFMpegWriter

writer = FFMpegWriter(fps=60, bitrate=8000)
mpl.rcParams['animation.ffmpeg_path'] = r"C:\\Users\\aumva\\OneDrive\\Documents\\ffmpeg-8.0.1-full_build\\bin\\ffmpeg.exe"

pi = np.pi

LINESTYLE_MAP = {
    "Planet": {
        "linestyle": "-.",
        "linewidth": 1.8,
        "alpha": 0.8
    },
    "dwarf": {
        "linestyle": "--",
        "linewidth": 1,
        "alpha": 0.6
    },
    "comet": {
        "linestyle": ":",
        "linewidth": 1,
        "alpha": 0.7
    }
}

def gregorian_to_julian(day, month, year):
    a = (14 - month)//12
    y = year + 4800 - a
    m = month + 12*a - 3
    JD = day + ((153*m+2)//5) + 365*y + y//4 - y//100 + y//400 - 32045
    return float(JD)

def days_from_date(day, month, year):
    return gregorian_to_julian(day, month, year) - 2451545.0

def julian_to_gregorian(JD):
    JD += 0.5
    Z = int(JD); F = JD - Z
    alpha = int((Z - 1867216.25)/36524.25) if Z >= 2299161 else 0
    A = Z + 1 + alpha - alpha//4 if Z >= 2299161 else Z
    B = A + 1524
    C = int((B - 122.1)/365.25)
    D = int(365.25*C)
    E = int((B - D)/30.6001)
    day = B - D - int(30.6001*E) + F
    month = E - 1 if E < 14 else E - 13
    year = C - 4716 if month > 2 else C - 4715
    return int(day), int(month), int(year)

# INPUTS

DD, MM, YY = 28, 1, 2061
cam_angle_deg = 30
FRAME_MODE = "b"
DYNAMIC_VIEW = "y0"
YEARS_TO_ANIMATE = 1
total_frames = 3600
FRAME_ROT = np.pi

start_days = days_from_date(DD, MM, YY)
TOTAL_DAYS_TO_ANIMATE = YEARS_TO_ANIMATE * 365.256363
dt = TOTAL_DAYS_TO_ANIMATE / total_frames

# CONSTANTS

c = 695700
s= 1
R_j = 1.0
M_sun = 1.989e30

# MASSES


M = {
    "Mercury":3.3011e23, "Venus":4.8675e24, "Earth":5.9722e24,
    "Mars":6.4171e23, "Jupiter":1.8982e27, "Saturn":5.6834e26,
    "Uranus":8.6810e25, "Neptune":1.02409e26,
    "Halley's Comet":2.2e14,"Ceres":9.393e20,"Vesta":2.5907e20,
    "Pallas": 2.04e20, "Pluto":1.3025e22, "Orcus":5.478e20, "Haumea":3.952e21, 
    "Makemake":2.69e21,"Eris":1.646e22, "Gonggong":1.75e21, "Sedna":2e21,
    "2012 VP113":4e19, "Leleākūhonua":2e19, "2023 KQ14": 7.5e18, 
    "2004 VN112": 2e19, "Planet Nine (Theorised)": 0,"(99942) Apophis": 6.1e10
    
}


# ORBITAL ELEMENTS 

moons = [
    {"name":"Mercury","a":57.91e6/c,"e":0.205,"R":2439*s/c,
     "color":"lightgray","omega":np.deg2rad(29.124),
     "Omega":np.deg2rad(48.331),"i":np.deg2rad(7.005),
     "L0":np.deg2rad(252.250905),"Ldot":np.deg2rad(4.09233445),
     "precession":np.deg2rad(1.5556),"epoch":2451545.0, "type":"Planet"},

    {"name":"Venus","a":108.21e6/c,"e":0.006772,"R":6051*s/c,
     "color":"peachpuff","omega":np.deg2rad(54.884),
     "Omega":np.deg2rad(76.680),"i":np.deg2rad(3.394),
     "L0":np.deg2rad(181.979801),"Ldot":np.deg2rad(1.60213034),
     "precession":np.deg2rad(1599.74/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Earth","a":149.59e6/c,"e":0.0167086,"R":6371*s/c,
     "color":"deepskyblue","omega":np.deg2rad(114.207),
     "Omega":np.deg2rad(-11.26064),"i":np.deg2rad(0.0),
     "L0":np.deg2rad(100.466456),"Ldot":np.deg2rad(0.98564736),
     "precession":np.deg2rad(1162.90/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Mars","a":227.93e6/c,"e":0.0934,"R":3389.5*s/c,
     "color":"red","omega":np.deg2rad(286.5),
     "Omega":np.deg2rad(49.558),"i":np.deg2rad(1.850),
     "L0":np.deg2rad(355.433000),"Ldot":np.deg2rad(0.52403900),
     "precession":np.deg2rad(1560.78/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Jupiter","a":778.479e6/c,"e":0.0489,"R":69911*s/c,
     "color":"tan","omega":np.deg2rad(273.867),
     "Omega":np.deg2rad(100.464),"i":np.deg2rad(1.303),
     "L0":np.deg2rad(34.351484),"Ldot":np.deg2rad(0.08308677),
     "precession":np.deg2rad(303.4/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Saturn","a":1433.53e6/c,"e":0.0565,"R":58232*s/c,
     "color":"goldenrod","omega":np.deg2rad(339.392),
     "Omega":np.deg2rad(113.665),"i":np.deg2rad(2.485),
     "L0":np.deg2rad(50.077471),"Ldot":np.deg2rad(0.03345965),
     "precession":np.deg2rad(150.9/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Uranus","a":2.870972e9/c,"e":0.04717,"R":25362*s/c,
     "color":"cyan","omega":np.deg2rad(96.9988),
     "Omega":np.deg2rad(74.006),"i":np.deg2rad(0.773),
     "L0":np.deg2rad(314.055005),"Ldot":np.deg2rad(0.01173129),
     "precession":np.deg2rad(75.25/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Neptune","a":4.50e9/c,"e":0.008678,"R":24622*s/c,
     "color":"cornflowerblue","omega":np.deg2rad(273.187),
     "Omega":np.deg2rad(131.784),"i":np.deg2rad(1.770),
     "L0":np.deg2rad(304.348665),"Ldot":np.deg2rad(0.00598103),
     "precession":np.deg2rad(38.13/3600),"epoch":2451545.0, "type":"Planet"},

    {"name":"Halley's Comet","a":17.737*149.59e6/c,"e":0.96658,"period":74.7*365.256363,"R":11*s/c,
     "color":"Mediumspringgreen","omega":np.deg2rad(112.05),
     "Omega":np.deg2rad(59.396),"i":np.deg2rad(161.96),
     "L0":np.deg2rad(171.51923),"Ldot":np.deg2rad(0.013087),
     "precession":np.deg2rad(0.002/3600),"epoch":2474040.5, "type":"comet"},

    {"name":"Ceres","a":413.7e6/c,"e":0.0785,"R":469.7*s/c,
     "color":"Magenta","omega":np.deg2rad(73.5977),
     "Omega":np.deg2rad(80.3055),"i":np.deg2rad(10.594),
     "L0":np.deg2rad(85.3032),"Ldot":np.deg2rad(0.214033),
     "precession":np.deg2rad(5400/3600),"epoch":2459600.5, "type":"dwarf"},

    {"name":"Vesta","a":353.4e6/c,"e":0.0894,"R":262.7*s/c,
     "color":"lightpink","omega":np.deg2rad(151.66),
     "Omega":np.deg2rad(103.71),"i":np.deg2rad(7.1422),
     "L0":np.deg2rad(64.4501),"Ldot":np.deg2rad(0.271834),
     "precession":np.deg2rad(3687/3600),"epoch":2453300.5, "type":"comet"},
    
    {"name":"Pallas","a":414.01e6/c,"e":0.2302,"R":256*s/c,
     "color":"deeppink","omega":np.deg2rad(310.9),
     "Omega":np.deg2rad(172.9),"i":np.deg2rad(34.93),
     "L0":np.deg2rad(164.4),"Ldot":np.deg2rad(0.21378),
     "precession":np.deg2rad(1700/3600),"epoch":2453300.5, "type":"comet"},
    
    {"name":"Pluto","a":5.90638e9/c,"e":0.2488,"R":1188*s/c,
     "color":"gainsboro","omega":np.deg2rad(113.834),
     "Omega":np.deg2rad(110.299),"i":np.deg2rad(17.16),
     "L0":np.deg2rad(238.5968),"Ldot":np.deg2rad(0.0039760),
     "precession":np.deg2rad(140/3600),"epoch":2451545.0, "type":"dwarf"},
    
    {"name":"Orcus","a":5.8606e9/c,"e":0.22701,"R":450*s/c,
     "color":"slategrey","omega":np.deg2rad(72.310),
     "Omega":np.deg2rad(268.799),"i":np.deg2rad(20.592),
     "L0":np.deg2rad(162.844),"Ldot":np.deg2rad(0.0039760),
     "precession":np.deg2rad(0/3600),"epoch":2459000.5, "type":"dwarf"},
    
    {"name":"Haumea","a":6.45e9/c,"e":0.19642,"period":103410,"R":780*s/c,
     "color":"greenyellow","omega":np.deg2rad(239.041),
     "Omega":np.deg2rad(121.8),"i":np.deg2rad(28.2137),
     "L0":np.deg2rad(219.413),"Ldot":np.deg2rad(0.00348),
     "precession":np.deg2rad(75/3600),"epoch":2459200.5, "type":"dwarf"},

    {"name":"Makemake","a":6.806e9/c,"e":0.1604,"period":112022,"R":715*s/c,
     "color":"yellow","omega":np.deg2rad(296.065),
     "Omega":np.deg2rad(79.441),"i":np.deg2rad(29.002),
     "L0":np.deg2rad(185.867),"Ldot":np.deg2rad(0.0032),
     "precession":np.deg2rad(75/3600),"epoch":2461000.5, "type":"dwarf"},

    {"name":"Eris","a":10.171e9/c,"e":0.4370,"R":1163*s/c,
     "color":"orange","omega":np.deg2rad(150.73),
     "Omega":np.deg2rad(36.02),"i":np.deg2rad(43.86),
     "L0":np.deg2rad(38.1310),"Ldot":np.deg2rad(0.0017661),
     "precession":np.deg2rad(59.81/3600),"epoch":2461000.5, "type":"dwarf"},
    
    {"name":"Gonggong","a":10.0073e9/c,"e":0.50318,"R":615*s/c,
     "color":"tomato","omega":np.deg2rad(206.6442),
     "Omega":np.deg2rad(336.8401),"i":np.deg2rad(30.8664),
     "L0":np.deg2rad(294.8683),"Ldot":np.deg2rad(0.001802),
     "precession":np.deg2rad(0/3600),"epoch":2461000.5, "type":"dwarf"},

    {"name":"Sedna","a":76e9/c,"e":0.8496,"R":995*s/c,
     "color":"hotpink","omega":np.deg2rad(311.352),
     "Omega":np.deg2rad(144.5),"i":np.deg2rad(11.93),
     "L0":np.deg2rad(93.717),"Ldot":np.deg2rad(0.0001),
     "precession":np.deg2rad(10/3600),"epoch":2458900.5, "type":"dwarf"},
    
    {"name":"2012 VP113","a":262.3*149.59e6/c,"e":0.6931,"R":225*s/c,
     "color":"violet","omega":np.deg2rad(293.90),
     "Omega":np.deg2rad(90.8),"i":np.deg2rad(24.0563),
     "L0":np.deg2rad(48.75),"Ldot":np.deg2rad(0.000232),
     "precession":np.deg2rad(10/3600),"epoch":2460800.5, "type":"dwarf"},
    
    {"name":"Leleākūhonua","a":1193*149.59e6/c,"e":0.94572,"R":110*s/c,
     "color":"mediumslateblue","omega":np.deg2rad(118.236),
     "Omega":np.deg2rad(300.989),"i":np.deg2rad(11.671),
     "L0":np.deg2rad(58.74),"Ldot":np.deg2rad(0.00002391),
     "precession":np.deg2rad(10/3600),"epoch":2460000.5, "type":"dwarf"},
    
    {"name":"2023 KQ14","a":252*149.59e6/c,"e":0.7385,"R":120*s/c,
     "color":"thistle","omega":np.deg2rad(198.74),
     "Omega":np.deg2rad(72.10),"i":np.deg2rad(10.98),
     "L0":np.deg2rad(267.4),"Ldot":np.deg2rad(0.00024652),
     "precession":np.deg2rad(10/3600),"epoch":2460800.5, "type":"dwarf"},
    
    {"name":"2004 VN112","a":608*149.59e6/c,"e":0.8579,"R":150*s/c,
     "color":"thistle","omega":np.deg2rad(326.72),
     "Omega":np.deg2rad(65.996),"i":np.deg2rad(25.572),
     "L0":np.deg2rad(33.398),"Ldot":np.deg2rad(0.000167),
     "precession":np.deg2rad(10/3600),"epoch":2459200.5, "type":"dwarf"},
    
    {"name":"(99942) Apophis","a":0.9224*149.59e6/c,"e":0.1911,"R":10*s/c,
     "color":"orange","omega":np.deg2rad(126.7),
     "Omega":np.deg2rad(203.9),"i":np.deg2rad(3.341),
     "L0":np.deg2rad(60.88),"Ldot":np.deg2rad(1.112),
     "precession":np.deg2rad(0/3600),"epoch":2460800.5, "type":"comet"},
    
    
]


# ORBIT CALCULATION 

def rotate_2d(x, y, theta):
    c, s = np.cos(theta), np.sin(theta)
    return c*x - s*y, s*x + c*y

def omega_at_time(m, days):
    """Argument of perihelion ω(t)"""
    JD = 2451545.0 + days
    centuries = (JD - m["epoch"]) / 36525.6363
    return m["omega"] + m["precession"] * centuries

def varpi_at_time(m, days):
    """Longitude of perihelion ϖ = Ω + ω"""
    return m["Omega"] + omega_at_time(m, days)

def mean_anomaly(days, m):
    JD = 2451545.0 + days
    L = m["L0"] + m["Ldot"] * (JD - m["epoch"])
    return np.mod(L - varpi_at_time(m, days), 2*np.pi)
    
def solve_kepler_elliptic(M, e, tol=1e-12):
    M = np.mod(M, 2*np.pi)
    E = M if e < 0.8 else np.pi

    for _ in range(50):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E

def compute_orbit_3d(days, m):
    M = mean_anomaly(days, m)
    e, a = m["e"], m["a"]
   
    if e < 1.0:
        E = M
        for _ in range(10):  
            E -= (E - e*np.sin(E) - M) / (1 - e*np.cos(E))

        nu = 2*np.arctan2(
            np.sqrt(1 + e) * np.sin(E / 2),
            np.sqrt(1 - e) * np.cos(E / 2)
        )

        r = a * (1 - e*np.cos(E))

    else:
        H = np.arcsinh(M / e)  
        for _ in range(10):
            H -= (e*np.sinh(H) - H - M) / (e*np.cosh(H) - 1)

        nu = 2*np.arctan2(
            np.sqrt(e + 1) * np.sinh(H / 2),
            np.sqrt(e - 1)
        )

        r = a * (e*np.cosh(H) - 1)

    x_op = r * np.cos(nu)
    y_op = r * np.sin(nu)

    ω = omega_at_time(m, days)
    cω, sω = np.cos(ω), np.sin(ω)

    x1 =  cω*x_op - sω*y_op
    y1 =  sω*x_op + cω*y_op

    i = m["i"]
    ci, si = np.cos(i), np.sin(i)

    x2 = x1
    y2 =  ci*y1
    z2 =  si*y1

    Ω = m["Omega"]
    cO, sO = np.cos(Ω), np.sin(Ω)

    X = cO*x2 - sO*y2
    Y = sO*x2 + cO*y2
    Z = z2

    return -X, -Y, Z

def project_2d(X, Y, Z, view_angle):
    # camera tilt
    ca, sa = np.cos(view_angle), np.sin(view_angle)
    Yp = ca*Y - sa*Z

    x2d, y2d = X, Yp

    
    return rotate_2d(x2d, y2d, FRAME_ROT)

def body_velocity_3d(days, m, dt=0.001):
    """
    Returns barycentric velocity of body in km/s
    days : simulation time (days)
    dt   : timestep in days (default = 1 day)
    """

    # positions in km
    x1, y1, z1 = compute_orbit_3d(days, m)
    x2, y2, z2 = compute_orbit_3d(days + dt, m)

    # km/day
    vx = (x2 - x1) / (dt)
    vy = (y2 - y1) / (dt)
    vz = (z2 - z1) / (dt)

    # convert km/day → km/s
    km_per_day_to_km_per_s = c / 86400.0
    vx *= km_per_day_to_km_per_s
    vy *= km_per_day_to_km_per_s
    vz *= km_per_day_to_km_per_s

    return vx, vy, vz

def draw_ecliptic(ax, view_angle,label="Ecliptic"):

    
    xmin, xmax = ax.get_xlim()
    L = max(abs(xmin), abs(xmax))
    X = np.array([-L, L])
    Y = np.array([0.0, 0.0])
    Z = np.array([0.0, 0.0])

    xe, ye = project_2d(X, Y, Z, view_angle)
    ax.plot(
            xe, ye,
            color="white",
            lw=1,
            alpha=0.5,
            zorder=0  
                     )

    ax.arrow(
        -xe[0], ye[0],
        -(xe[1] - xe[0]), -(ye[1] - ye[0]),
        head_width=0.04 * L,
        head_length=0.06 * L,
        color="white",
        linewidth=1,
        length_includes_head=True,
        alpha=0.6
    )

    xl =  xe[1]+0.28*L
    yl =  ye[1]-0.05*L
    ax.text(
        -xl, yl,
        "Vernal Equinox",
        color="white",
        ha="left",
        va="center",
        fontsize=10
    )

    X = np.array([0.0, 0.0])
    Y = np.array([-L, L])
    Z = np.array([0.0, 0.0])

    xe, ye = project_2d(X, Y, Z, view_angle)

    ax.plot(
        xe, ye,
        color="white",
        linestyle= "--",
        linewidth=1.5,
        alpha=0.6
    )

    tx, ty = project_2d(L * 0.85, -L*0.05, 0.0, view_angle)

    ax.text(
    tx, ty,
        label,
        color="white",
        fontsize=10,
        ha="right",
        va="center",
        alpha=1
    )

def peri_apo_3d(days, m):
    
    """
    Returns 3D Cartesian coordinates of
    perihelion (q) and aphelion (Q)
    """

    a = m["a"]
    e = m["e"]


    r_peri = a * (1 - e)
    r_apo  = a * (1 + e)


    xq, yq = r_peri, 0.0
    xQ, yQ = -r_apo, 0.0


    ω = omega_at_time(m, days)
    cω, sω = np.cos(ω), np.sin(ω)

    xq1 =  cω*xq - sω*yq
    yq1 =  sω*xq + cω*yq

    xQ1 =  cω*xQ - sω*yQ
    yQ1 =  sω*xQ + cω*yQ

    i = m["i"]
    ci, si = np.cos(i), np.sin(i)

    xq2 = xq1
    yq2 = ci*yq1
    zq2 = si*yq1

    xQ2 = xQ1
    yQ2 = ci*yQ1
    zQ2 = si*yQ1

    Ω = m["Omega"]
    cO, sO = np.cos(Ω), np.sin(Ω)

    Xq =  cO*xq2 - sO*yq2
    Yq =  sO*xq2 + cO*yq2
    Zq =  zq2

    XQ =  cO*xQ2 - sO*yQ2
    YQ =  sO*xQ2 + cO*yQ2
    ZQ =  zQ2

    return (-Xq, -Yq, Zq), (-XQ, -YQ, ZQ)

def node_positions_3d(m, days, npts=400):
    
    nu = np.linspace(0, 2*np.pi, npts)
    a, e = m["a"], m["e"]

    r = a * (1 - e*e) / (1 + e*np.cos(nu))

    x_op = r * np.cos(nu)
    y_op = r * np.sin(nu)
    z_op = np.zeros_like(x_op)

    ω = omega_at_time(m, days)
    cω, sω = np.cos(ω), np.sin(ω)

    x1 =  cω*x_op - sω*y_op
    y1 =  sω*x_op + cω*y_op

    i = m["i"]
    ci, si = np.cos(i), np.sin(i)

    x2 = x1
    y2 =  ci*y1
    z2 =  si*y1

    Ω = m["Omega"]
    cO, sO = np.cos(Ω), np.sin(Ω)

    X =  cO*x2 - sO*y2
    Y =  sO*x2 + cO*y2
    Z =  z2

    asc = None
    desc = None

    for k in range(len(Z)-1):
        if Z[k] == 0:
            continue

        if Z[k] < 0 and Z[k+1] > 0:
            t = -Z[k] / (Z[k+1] - Z[k])
            asc = (
                X[k] + t*(X[k+1]-X[k]),
                Y[k] + t*(Y[k+1]-Y[k]),
                0.0
            )

        if Z[k] > 0 and Z[k+1] < 0:
            t = -Z[k] / (Z[k+1] - Z[k])
            desc = (
                X[k] + t*(X[k+1]-X[k]),
                Y[k] + t*(Y[k+1]-Y[k]),
                0.0
            )

    return asc, desc

def depth_style(Z, base_alpha, base_lw):
    Z = np.asarray(Z, dtype=float)
    mask_above = Z > 0
    mask_below = Z < 0

    return {
        "above": dict(alpha=base_alpha, linewidth=base_lw),
        "below": dict(alpha=base_alpha * 0.4, linewidth=base_lw),
        "mask_above": mask_above,
        "mask_below": mask_below
    }

def barycentric_vectors(days):
    vectors = {}
    total_mass = M_sun + sum(M[m["name"]] for m in moons)
    GIANTS = ["Jupiter", "Saturn", "Uranus", "Neptune"]

    for m in moons:
        if m["name"] not in GIANTS:
            continue

        X, Y, Z = compute_orbit_3d(days, m)
        mass = M[m["name"]]

        vectors[m["name"]] = (
            mass * X / total_mass,
            mass * Y / total_mass,
            mass * Z / total_mass
        )

    return vectors


# BARYCENTER ANIMATION 

fig_zoom, ax_zoom = plt.subplots(figsize=(16,16))

moving_x = []
moving_y = []
moving_z = []
heli_X = []
heli_Y = []
heli_Z = []
bary_dist = []
bary_time = []

def gravity_barycenter_3d(days):
    """Compute full 3D barycenter (Sun at origin)"""
    total_mass = M_sun + sum(M[m["name"]] for m in moons)
    
    BX = BY = BZ = 0.0
    for m in moons:
        X, Y, Z = compute_orbit_3d(days, m)
        mass = M[m["name"]]
        BX += mass * X
        BY += mass * Y
        BZ += mass * Z

    return BX / total_mass, BY / total_mass, BZ / total_mass

def draw_zoomb(ax, days):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    zoom = 3
    try:
        
        BX, BY, BZ = gravity_barycenter_3d(days)
        true_dist = np.sqrt(BX**2 + BY**2 + BZ**2)

        if FRAME_MODE == "h":
           Xb, Yb, Zb = BX, BY, BZ
           Xs, Ys, Zs = 0.0, 0.0, 0.0

        elif FRAME_MODE == "b":
             Xb, Yb, Zb = 0.0, 0.0, 0.0
             Xs, Ys, Zs = -BX, -BY, -BZ

        if FRAME_MODE == "h":
           moving_x.append(Xb)
           moving_y.append(Yb)
           moving_z.append(Zb)
           
        elif FRAME_MODE == "b":
             moving_x.append(Xs)
             moving_y.append(Ys)
             moving_z.append(Zs)

        bary_dist.append(true_dist)
        bary_time.append(days)

        mx,my,mz = np.array(moving_x), np.array(moving_y), np.array(moving_z)
        mx_above, my_above, mz_above = mx.copy(), my.copy(), mz.copy()
        mx_above[mz < 0], my_above[mz < 0], mz_above[mz < 0] = np.nan, np.nan, np.nan
        mx_below, my_below, mz_below = mx.copy(), my.copy(), mz.copy()
        mx_below[mz >= 0], my_below[mz >= 0], mz_below[mz >= 0] = np.nan, np.nan, np.nan
        pxa, pya = project_2d(mx_above, my_above, mz_above, view_angle)
        pxb, pyb = project_2d(mx_below, my_below, mz_below, view_angle)
        ax.plot(pxa, pya, "-.", color="red", linewidth=2, alpha=0.9)
        ax.plot(pxb, pyb, "-.", color="red", linewidth=2, alpha=0.4)
        hx, hy = project_2d(mx[-1], my[-1], mz[-1], view_angle)
        bx, by = project_2d(Xb, Yb, Zb, view_angle)

        if FRAME_MODE == "h":          
            ax.add_patch(plt.Circle((0, 0), R_j, color="yellow"))
            ax.text(0, -(R_j+0.2), "Sun", color="white", ha="center")
            ax.plot(0, 0, "x", color="black", markersize=10, mew=2)
            ax.text(0, -0.3, "Heliocenter", ha="center", color="black")
            
            alpha_b = 1.0 if Zb >= 0 else 0.4
            ax.plot(bx, by, "x", color="red", markersize=10, mew=2, alpha=alpha_b)
            ax.text(bx, by - 0.3, "Barycenter",
                    ha="center", color="red")

        elif FRAME_MODE == "b":
            ax.plot(0, 0, "x", color="red", markersize=15, mew=2)
            ax.text(0, -0.3, "Barycenter", ha="center", color="red")          

            vt = 0.01
            BX1, BY1, BZ1 = gravity_barycenter_3d(days)
            BX2, BY2, BZ2 = gravity_barycenter_3d(days + vt)
            BVX = 1000*c*(BX2 - BX1) / (vt*86400)
            BVY = 1000*c*(BY2 - BY1) / (vt*86400)
            BVZ = 1000*c*(BZ2 - BZ1) / (vt*86400)
            VSX, VSY, VSZ = -BVX, -BVY, -BVZ
            bvmag = np.sqrt(VSX**2 + VSY**2 + VSZ**2)
            vx2d, vy2d = project_2d(VSX, VSY, VSZ, view_angle)

            vectors = barycentric_vectors(days)
            rpv_contrib = {"Jupiter": 0.0, "Saturn": 0.0, "Uranus": 0.0,"Neptune": 0.0}

            for name, (VX, VY, VZ) in vectors.items():
                BX, BY, BZ = gravity_barycenter_3d(days)
                bx, by = project_2d(BX, BY, BZ, view_angle)
                vx, vy = project_2d(VX, VY, VZ, view_angle)
                mag = np.sqrt(VX**2 + VY**2 + VZ**2)
                rmag = np.sqrt(BX**2 + BY**2 + BZ**2)
                rpv_contrib[name] += mag
                alpha_s1 = 1.0 if mz[-1] >= 0 else 0.2
                alpha_s2 = 1.0 if BZ >= 0 else 0.2
                alpha_s3 = 1.0 if VSZ >= 0 else 0.1
                alpha_s  = 1.0 if VZ >= 0 else 0.4

                PLANET_COLOR = {m["name"]:m["color"] for m in moons} 

                ax.arrow(0, 0, vx, vy, length_includes_head=True,
                         head_width=0.06, head_length=0.12,
                         color=PLANET_COLOR[name], linewidth=1.5,
                         alpha=alpha_s, zorder=6)           
                ax.arrow(0, 0, bx, by, color="lightgray", linewidth=2,
                         head_width=0.06, head_length=0.12,
                         alpha=alpha_s2, zorder=7)
                ax.arrow(hx, hy, 0.04*vx2d, 0.04*vy2d,
                         length_includes_head=True,
                         head_width=0.03, head_length=0.06,
                         color="deepskyblue", linewidth=2,
                         alpha=alpha_s3, zorder=8)

                ax.text(bx*1.2, by*1.2, "Resultant Pull Vector",
                        color="lightgray", fontsize=10,
                        ha="center", va="center", alpha=alpha_s2)

                ax.text(vx*1.1, vy*1.1, name,
                        color=PLANET_COLOR[name],
                        fontsize=9, ha="center",
                        va="center", alpha=alpha_s)

                ax.add_patch(plt.Circle((hx, hy), R_j,
                                         color="yellow", alpha=alpha_s1))
                ax.text(hx, hy - (R_j+0.2), "Sun",
                        color="white", ha="center", alpha=alpha_s1)
                ax.plot(hx, hy, "x", color="black",
                        markersize=15, mew=2, alpha=alpha_s1)
                ax.text(hx, hy - 0.3, "Heliocenter",
                        ha="center", color="black", alpha=alpha_s1)

                
                
        JD_now = 2451545.0 + days
        d0, m0, y0 = julian_to_gregorian(JD_now)
        
        if YEARS_TO_ANIMATE == 0 or TOTAL_DAYS_TO_ANIMATE == 0:
               ax.set_title(
                f"Barycentric Deflection of the Sun at {d0:02d}/{m0:02d}/{y0}\n"
                f"Viewing Angle = {cam_angle_deg:.0f} Degrees from the Ecliptic\n"
                f"Current distance = {true_dist:.3f} Solar Radii\n"
                f"Dimming = Below the Ecliptic\n"
                f"Brightening = Above the Ecliptic",
                color="black"
                )
               ax.text(
                0.98, 0.98,               
                f"Heliocentric Deflection = {mag:.3f} Solar Radii\n"
                f"Barycentric Velocity = {bvmag:.2f} m/s",
                transform=ax.transAxes,
                ha="right", va="top",
                fontsize=10, color="cyan",
                bbox=dict(facecolor="black", edgecolor="cyan", alpha=0.7)
               )
               ax.text(
                       0.02, 0.98,
                       f"Resultant Pull Vector Deflection = {mag:.3f} Solar Radii\n"
                       f"Jupiter Deflection= {-rpv_contrib['Jupiter']:.3f} Solar Radii\n"
                       f"Saturn Deflection= {-rpv_contrib['Saturn']:.3f} Solar Radii\n"
                       f"Uranus Deflection= {-rpv_contrib['Uranus']:.3f} Solar Radii\n"
                       f"Neptune Deflection= {-rpv_contrib['Neptune']:.3f} Solar Radii",
                       transform=ax.transAxes,
                       ha="left", va="top",
                       fontsize=10, color="orange",
                       bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                         )

        else:
                 max_dist = max(bary_dist)
                 max_idx = int(np.argmax(bary_dist))
                 max_x = moving_x[max_idx]
                 max_y = moving_y[max_idx]
                 max_z = moving_z[max_idx]
                 mbx, mby = project_2d(max_x,max_y,max_z,view_angle)
                 max_days = bary_time[max_idx]
                 ax.plot(mbx, mby, "o", color="cyan", markersize=5)
                 d_max, m_max, y_max = julian_to_gregorian(2451545.0 + max_days)
                 ax.set_title(
                              f"Barycentric wobble of the Sun\n"
                              f"{DD}/{MM}/{YY} to {d0:02d}/{m0:02d}/{y0}\n"
                              f"Simulation time = {YEARS_TO_ANIMATE} Years\n"
                              f"Viewing Angle = {cam_angle_deg:.0f} Degrees from the Ecliptic\n"
                              f"Dimming = Below the Ecliptic\n"
                              f"Brightening = Above the Ecliptic",
                              color="black"
                 )
                 
            
                 ax.text(
                         0.98, 0.98,
                         f"Max Heliocentric Deflection = {max_dist:.3f} Solar Radii\n"
                         f"Max Heliocentric Deflection at = {d_max:02d}/{m_max:02d}/{y_max}\n"
                         f"Heliocentric Deflection = {rmag:.3f} Solar Radii\n"
                         f"Barycentric Velocity = {bvmag:.2f} m/s",
                         transform=ax.transAxes,
                         ha="right", va="top",
                         fontsize=10, color="cyan",
                         bbox=dict(facecolor="black", edgecolor="cyan", alpha=0.7)
                         )
                 ax.text(
                         0.02, 0.98,
                         f"Resultant Pull Vector Deflection= {-rmag:.3f} Solar Radii\n"
                         f"Jupiter Deflection = {-rpv_contrib['Jupiter']:.3f} Solar Radii\n"
                         f"Saturn Deflection = {-rpv_contrib['Saturn']:.3f} Solar Radii\n"
                         f"Uranus Deflection = {-rpv_contrib['Uranus']:.3f} Solar Radii\n"
                         f"Neptune Deflection = {-rpv_contrib['Neptune']:.3f} Solar Radii",
                         transform=ax.transAxes,
                         ha="left", va="top",
                         fontsize=10, color="orange",
                         bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                           )
            
        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(-zoom, zoom)
        ticks = np.linspace(-zoom, zoom, 9)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:.3f}" for t in ticks], color="black")
        ax.set_yticklabels([f"{t:.3f}" for t in ticks], color="black")
        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(-zoom, zoom)
        draw_ecliptic(ax, view_angle)
        ax.set_xlabel("Distance in Solar Radii (1 Solar Radii = 695,700 km)", color="black")
        ax.set_ylabel("Distance in Solar Radii (1 Solar Radii = 695,700 km)", color="black")
        ax.grid(True, linestyle=":", linewidth=0.8, color="white", alpha=0.5)
        ax.set_axisbelow(True)

    except Exception as exc:
        ax.set_facecolor("black")
        ax.text(
            0.5, 0.5,
            f"Error in draw_zoomb:\n{type(exc).__name__}: {str(exc)}",
            color="white",
            ha="center",
            va="center",
            transform=ax.transAxes
        )
        ax.set_xticks([])
        ax.set_yticks([])
        
def update_bary(frame):
    global cam_angle_deg, cam_angle, view_angle
    if DYNAMIC_VIEW == "y":
        frac = (frame / total_frames)**0.5
        frac = min(frac, 1.0)                 
        cam_angle_deg = 90 * (1 - frac)      
        cam_angle_deg = max(cam_angle_deg, 0) 
        cam_angle = np.deg2rad(cam_angle_deg)
        view_angle = np.pi/2 - cam_angle
    else:
         cam_angle = np.deg2rad(cam_angle_deg)
         view_angle = np.pi/2 - cam_angle
    
    draw_zoomb(ax_zoom, start_days + frame*dt)
ani_zoom = animation.FuncAnimation(fig_zoom, update_bary,
                                        frames=total_frames, interval=17)
plt.show()

# INNER PLANETS

fig_inner, ax_inner = plt.subplots(figsize=(24,24))

def draw_inner(ax, days):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    amax = 227.93e6/ c                                        #413.7e6
    inner = [m for m in moons if m["name"] in
             ["Mercury", "Venus", "Earth", "Mars", "Jupiter","Saturn",
              "Uranus","Neptune", "Halley's Comet","Ceres","Vesta","Pallas"]]

    for m in inner:
        X, Y, Z = compute_orbit_3d(days, m)
        xp, yp = project_2d(X, Y, Z, view_angle)
        ax.add_patch(plt.Circle((xp, yp), m["R"], color=m["color"]))
        ax.text(
            xp, yp,
            m["name"],
            color=m["color"],
            ha="center",
            va="center",
            fontsize=10,
            zorder=1,
            clip_on=True
        )
        
        nu_ell = np.linspace(0, 2*pi, 400)
        Xo, Yo, Zo = [], [], []
        for nu in nu_ell:
            r = m["a"] * (1 - m["e"]**2) / (1 + m["e"] * np.cos(nu))
            x = r * np.cos(nu)
            y = r * np.sin(nu)

            ω = omega_at_time(m, days)
            cω, sω = np.cos(ω), np.sin(ω)
            x1 =  cω*x - sω*y
            y1 =  sω*x + cω*y

            ci, si = np.cos(m["i"]), np.sin(m["i"])
            x2 = x1
            y2 = ci*y1
            z2 = si*y1

            Omega_t = m["Omega"]
            cO, sO = np.cos(Omega_t), np.sin(Omega_t)
            Xo.append(cO*x2 - sO*y2)
            Yo.append(sO*x2 + cO*y2)
            Zo.append(z2)
                
        xx, yy = project_2d(-np.array(Xo), -np.array(Yo), np.array(Zo), view_angle)
        Z = np.array(Zo)
        xx_above, yy_above= xx.copy(), yy.copy()
        xx_below, yy_below = xx.copy(), yy.copy()
        xx_above[Z < 0], yy_above[Z < 0]  = np.nan, np.nan
        xx_below[Z >= 0], yy_below[Z >= 0] = np.nan, np.nan
        ls = LINESTYLE_MAP.get(m["type"])
        style = depth_style(Zo, ls["alpha"], ls["linewidth"])
        ax.plot(xx_above, yy_above,
        linestyle=ls["linestyle"],
        color=m["color"],
        alpha=style["above"]["alpha"],
        linewidth=style["above"]["linewidth"])
        ax.plot(xx_below, yy_below,
        linestyle=ls["linestyle"],
        color=m["color"],
        alpha=style["below"]["alpha"],
        linewidth=style["below"]["linewidth"])

        if m["e"] < 1:   
          (Xp, Yp, Zp), (Xa, Ya, Za) = peri_apo_3d(days, m)
          px, py = project_2d(Xp, Yp, Zp, view_angle)
          ax_, ay_ = project_2d(Xa, Ya, Za, view_angle)
          ax.plot(px, py, "o", color=m["color"], markersize=2, zorder=5)
          ax.text(px, py - 10*0.5, "q",
            color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
          ax.plot(ax_, ay_, "o", color=m["color"], markersize=2, zorder=5)
          ax.text(ax_, ay_ - 10*0.5, "Q",
            color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
          asc, desc = node_positions_3d(m, days)

          if asc is not None:
             ax_, ay_ = project_2d(*asc, view_angle)
             ax.plot(-ax_, -ay_, "s", color=m["color"], markersize=2)
             ax.text(-ax_, -(ay_ + 10*0.5), "☊", color=m["color"],
                     fontsize=8, ha="center", va="top", clip_on=True)

          if desc is not None:
             dx_, dy_ = project_2d(*desc, view_angle)
             ax.plot(-dx_, -dy_, "s", color=m["color"], markersize=2)
             ax.text(-dx_, -(dy_ + 10*0.5), "☋", color=m["color"],
            fontsize=8, ha="center", va="top", clip_on=True)
        
    #ax.plot(0, 0, "x", color="red", zorder=10)
    ax.add_patch(plt.Circle((0, 0), R_j, color="yellow", zorder=9))
    ax.text(0, -20, "Sun", color="yellow", ha="center")
    peri, apo = peri_apo_3d(days, m)
    
    ax.set_xlim(-1.01 * amax*(1+0.2302), 1.01 * amax*(1+0.2302))
    ax.set_ylim(-1.01 * amax*(1+0.2302), 1.01 * amax*(1+0.2302))
    draw_ecliptic(ax, view_angle)
    ticks = np.linspace(-1.01 * amax*(1+0.2302), 1.01 * amax*(1+0.2302), 9)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_yticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_xlabel("Distance in AU (1AU = 149,597,870 km)")
    ax.set_ylabel("Distance in AU (1AU = 149,597,870 km)")
    ax.grid(True, linestyle=":", linewidth=0.8, color="white", alpha=0.3)
    ax.set_axisbelow(True)
    velocities = {}       
    for m in inner:   
        vx_v, vy_v, vz_v = body_velocity_3d(days, m)
        speed = np.sqrt(vx_v**2 + vy_v**2 + vz_v**2)
        velocities[m["name"]] = speed

    lines = ["Orbital Velocities (km/s)"]
    for name, v in sorted(velocities.items(), key=lambda x: -x[1]):
        lines.append(f"{name:8s} : {v:5.2f}")
        vel_text = "\n".join(lines)

    ax.text(0.98, 0.98,vel_text,transform=ax.transAxes,ha="right",va="top",fontsize=10,
            color="cyan", bbox=dict(facecolor="black",edgecolor="cyan",alpha=0.5,zorder=20))
    JD = 2451545.0 + days
    d, mth, y = julian_to_gregorian(JD)

    if YEARS_TO_ANIMATE == 0:
        ax.set_title(
            f"Inner Solar System at {d:02d}/{mth:02d}/{y}\n"
            f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n",
            color="black"
        )
        
        ax.text(
                0.02, 0.98,   
                "Aphelion = Q\n"
                "Perihelion = q\n"
                "Ascending Node = ☊\n"
                "Descending Node = ☋\n"
                "Above Ecliptic = Bright\n"
                "Below Ecliptic = Dim"
                ,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10, color="orange",
                bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                  )
    else:
        ax.set_title(
            f"Inner Solar System animation from {DD}/{MM}/{YY} to {d:02d}/{mth:02d}/{y}\n"
            f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n"
            f"Simulation time = {YEARS_TO_ANIMATE:.0f} Year",
            color="black"
        )
        ax.text(
                0.02, 0.98,   
                "Aphelion = Q\n"
                "Perihelion = q\n"
                "Ascending Node = ☊\n"
                "Descending Node = ☋\n"
                "Above Ecliptic = Bright\n"
                "Below Ecliptic = Dim"
                ,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10, color="orange",
                bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                  )


def update_inner(frame):
    global cam_angle_deg, cam_angle, view_angle
    if DYNAMIC_VIEW == "y":
        frac = (frame / total_frames)**0.5
        frac = min(frac, 1.0)                 
        cam_angle_deg = 90 * (1 - frac)      
        cam_angle_deg = max(cam_angle_deg, 0) 
        cam_angle = np.deg2rad(cam_angle_deg)
        view_angle = np.pi/2 - cam_angle
    else:
         cam_angle = np.deg2rad(cam_angle_deg)
         view_angle = np.pi/2 - cam_angle
    
    draw_inner(ax_inner, start_days + frame * dt)
ani_inner = animation.FuncAnimation(
    fig_inner,
    update_inner,
    frames=total_frames,
    interval=17
)

plt.show()


# OUTER PLANETS


fig_outer, ax_outer = plt.subplots(figsize=(16,16))

def draw_outer(ax, days):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    
    outer = [m for m in moons if m["name"] in
             ["Ceres","Jupiter","Saturn","Uranus","Neptune",
              "Pluto", "Orcus","Haumea","Makemake", 
               "Eris", "Gonggong", "Sedna", "Halley's Comet"]]
    amax = 4.50e9/c
    for m in outer:
       X, Y, Z = compute_orbit_3d(days, m)
       xp, yp = project_2d(X, Y, Z, view_angle)
       ax.add_patch(plt.Circle((xp, yp), m["R"], color=m["color"]))
       ax.text(
           xp, yp,
           m["name"],
           color=m["color"],
           ha="center",
           va="center",
           fontsize=10,
           zorder=6,
           clip_on=True
       )

       nu_ell = np.linspace(0, 2*pi, 400)
       Xo, Yo, Zo = [], [], []

       for nu in nu_ell:
           r = m["a"] * (1 - m["e"]**2) / (1 + m["e"] * np.cos(nu))
           x = r * np.cos(nu)
           y = r * np.sin(nu)

           ω = omega_at_time(m, days)
           cω, sω = np.cos(ω), np.sin(ω)
           x1 =  cω*x - sω*y
           y1 =  sω*x + cω*y

           ci, si = np.cos(m["i"]), np.sin(m["i"])
           x2 = x1
           y2 = ci*y1
           z2 = si*y1

           Omega_t = m["Omega"]
           cO, sO = np.cos(Omega_t), np.sin(Omega_t)
           Xo.append(cO*x2 - sO*y2)
           Yo.append(sO*x2 + cO*y2)
           Zo.append(z2)

       xx, yy = project_2d(-np.array(Xo), -np.array(Yo), np.array(Zo), view_angle)
       Z = np.array(Zo)
       xx_above, yy_above= xx.copy(), yy.copy()
       xx_below, yy_below = xx.copy(), yy.copy()
       xx_above[Z < 0], yy_above[Z < 0]  = np.nan, np.nan
       xx_below[Z >= 0], yy_below[Z >= 0] = np.nan, np.nan
       ls = LINESTYLE_MAP.get(m["type"])
       style = depth_style(Zo, ls["alpha"], ls["linewidth"])
       ax.plot(xx_above, yy_above,
       linestyle=ls["linestyle"],
       color=m["color"],
       alpha=style["above"]["alpha"],
       linewidth=style["above"]["linewidth"])
       ax.plot(xx_below, yy_below,
       linestyle=ls["linestyle"],
       color=m["color"],
       alpha=style["below"]["alpha"],
       linewidth=style["below"]["linewidth"])
       if m["e"] < 1:   
         (Xp, Yp, Zp), (Xa, Ya, Za) = peri_apo_3d(days, m)
         px, py = project_2d(Xp, Yp, Zp, view_angle)
         ax_, ay_ = project_2d(Xa, Ya, Za, view_angle)
         ax.plot(px, py, "o", color=m["color"], markersize=2, zorder=5)
         ax.text(px, py - 10, "q",
           color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
         ax.plot(ax_, ay_, "o", color=m["color"], markersize=2, zorder=5)
         ax.text(ax_, ay_ - 10, "Q",
           color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
         asc, desc = node_positions_3d(m, days)

         if asc is not None:
            ax_, ay_ = project_2d(*asc, view_angle)
            ax.plot(-ax_, -ay_, "s", color=m["color"], markersize=2)
            ax.text(-ax_, -(ay_ + 10), "☊", color=m["color"],
                    fontsize=8, ha="center", va="top", clip_on=True)

         if desc is not None:
            dx_, dy_ = project_2d(*desc, view_angle)
            ax.plot(-dx_, -dy_, "s", color=m["color"], markersize=2)
            ax.text(-dx_, -(dy_ + 10), "☋", color=m["color"],
           fontsize=8, ha="center", va="top", clip_on=True)
    #ax.plot(0, 0, "x", color="red", zorder=10)
    ax.add_patch(plt.Circle((0, 0), R_j, color="yellow", zorder=9))
    ax.text(0, -20, "Sun", color="yellow", ha="center")

    ax.set_xlim(-1.1*amax,1.1*amax)
    ax.set_ylim(-1.1*amax,1.1*amax)
    draw_ecliptic(ax, view_angle)
    ticks = np.linspace(-1.1*amax,1.1*amax,9)
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_yticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_xlabel("Distance in AU (1AU = 149,597,870 km)")
    ax.set_ylabel("Distance in AU (1AU = 149,597,870 km)")
    ax.grid(True, linestyle=':', linewidth=0.8, color='white', alpha=0.3)
    ax.set_axisbelow(True) 
    velocities = {}       
    for m in outer:   
        vx_v, vy_v, vz_v = body_velocity_3d(days, m)
        speed = np.sqrt(vx_v**2 + vy_v**2 + vz_v**2)
        velocities[m["name"]] = speed

    lines = ["Orbital Velocities (km/s)"]
    for name, v in sorted(velocities.items(), key=lambda x: -x[1]):
        lines.append(f"{name:8s} : {v:5.2f}")
        vel_text = "\n".join(lines)

    ax.text(0.98, 0.98,vel_text,transform=ax.transAxes,ha="right",va="top",fontsize=10,
            color="cyan", bbox=dict(facecolor="black",edgecolor="cyan",alpha=0.5,zorder=20))
    JD = 2451545.0 + days
    d, m, y = julian_to_gregorian(JD)
    
    if YEARS_TO_ANIMATE==0:
       ax.set_title(f"Outer Solar System at {d:02d}/{m:02d}/{y}\n"
                    f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n",
                    color="black")
       
       ax.text(
               0.02, 0.98,   
               "Aphelion = Q\n"
               "Perihelion = q\n"
               "Ascending Node = ☊\n"
               "Descending Node = ☋\n"
               "Above Ecliptic = Bright\n"
               "Below Ecliptic = Dim"
               ,
               transform=ax.transAxes,
               ha="left", va="top",
               fontsize=10, color="orange",
               bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                 )
    else:
        ax.set_title(f"Outer Solar System animation from {DD}/{MM}/{YY} to {d:02d}/{m:02d}/{y}\n"
                     f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n"
                     f"Simulation time = {YEARS_TO_ANIMATE:.0f} Years",
                     color="black")
        
        ax.text(
                0.02, 0.98,   
                "Aphelion = Q\n"
                "Perihelion = q\n"
                "Ascending Node = ☊\n"
                "Descending Node = ☋\n"
                "Above Ecliptic = Bright\n"
                "Below Ecliptic = Dim"
                ,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10, color="orange",
                bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                  )


def update_outer(frame):
    global cam_angle_deg, cam_angle, view_angle
    if DYNAMIC_VIEW == "y":
        frac = (frame / total_frames)**0.5
        frac = min(frac, 1.0)                 
        cam_angle_deg = 90 * (1 - frac)       
        cam_angle_deg = max(cam_angle_deg, 0) 
        cam_angle = np.deg2rad(cam_angle_deg)
        view_angle = np.pi/2 - cam_angle
    else:
         cam_angle = np.deg2rad(cam_angle_deg)
         view_angle = np.pi/2 - cam_angle
    
    draw_outer(ax_outer, start_days + frame*dt)
ani_outer = animation.FuncAnimation(fig_outer, update_outer,
                                        frames=total_frames, interval=17)
plt.show()


# KUIPER BELT 

fig_kb, ax_kb = plt.subplots(figsize=(16,16))

def draw_kb(ax, days):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_aspect("equal")
    amax = 10.171e9/c
    emax = 0.4370
    kb = [m for m in moons if m["name"] in
             ["Uranus", "Neptune","Pluto","Orcus","Haumea","Makemake", 
              "Eris", "Gonggong", "Sedna","Halley's Comet","2012 VP113"
              "Leleākūhonua","2023 KQ14","2004 VN112"]] 

    for m in kb:
        X, Y, Z = compute_orbit_3d(days, m)
        xp, yp = project_2d(X, Y, Z, view_angle)
        ax.add_patch(plt.Circle((xp, yp), m["R"], color=m["color"]))
        ax.text(
            xp, yp,
            m["name"],
            color=m["color"],
            ha="center",
            va="center",
            fontsize=10,
            zorder=1,
            clip_on=True
        )

        nu_ell = np.linspace(0, 2*pi, 400)
        Xo, Yo, Zo = [], [], []

        for nu in nu_ell:
            r = m["a"] * (1 - m["e"]**2) / (1 + m["e"] * np.cos(nu))
            x = r * np.cos(nu)
            y = r * np.sin(nu)

            ω = omega_at_time(m, days)
            cω, sω = np.cos(ω), np.sin(ω)
            x1 =  cω*x - sω*y
            y1 =  sω*x + cω*y

            ci, si = np.cos(m["i"]), np.sin(m["i"])
            x2 = x1
            y2 = ci*y1
            z2 = si*y1

            Omega_t = m["Omega"]
            cO, sO = np.cos(Omega_t), np.sin(Omega_t)
            Xo.append(cO*x2 - sO*y2)
            Yo.append(sO*x2 + cO*y2)
            Zo.append(z2)

        xx, yy = project_2d(-np.array(Xo), -np.array(Yo), np.array(Zo), view_angle)
        Z = np.array(Zo)
        xx_above, yy_above= xx.copy(), yy.copy()
        xx_below, yy_below = xx.copy(), yy.copy()
        xx_above[Z < 0], yy_above[Z < 0]  = np.nan, np.nan
        xx_below[Z >= 0], yy_below[Z >= 0] = np.nan, np.nan
        ls = LINESTYLE_MAP.get(m["type"])
        style = depth_style(Zo, ls["alpha"], ls["linewidth"])
        ax.plot(xx_above, yy_above,
        linestyle=ls["linestyle"],
        color=m["color"],
        alpha=style["above"]["alpha"],
        linewidth=style["above"]["linewidth"])
        ax.plot(xx_below, yy_below,
        linestyle=ls["linestyle"],
        color=m["color"],
        alpha=style["below"]["alpha"],
        linewidth=style["below"]["linewidth"])
        if m["e"] < 1:   
          (Xp, Yp, Zp), (Xa, Ya, Za) = peri_apo_3d(days, m)
          px, py = project_2d(Xp, Yp, Zp, view_angle)
          ax_, ay_ = project_2d(Xa, Ya, Za, view_angle)
          ax.plot(px, py, "o", color=m["color"], markersize=2, zorder=5)
          ax.text(px, py - 10, "q",
            color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
          ax.plot(ax_, ay_, "o", color=m["color"], markersize=2, zorder=5)
          ax.text(ax_, ay_ - 10, "Q",
            color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
          asc, desc = node_positions_3d(m, days)

          if asc is not None:
             ax_, ay_ = project_2d(*asc, view_angle)
             ax.plot(-ax_, -ay_, "s", color=m["color"], markersize=2)
             ax.text(-ax_, -(ay_ + 10), "☊", color=m["color"],
                     fontsize=8, ha="center", va="top", clip_on=True)

          if desc is not None:
             dx_, dy_ = project_2d(*desc, view_angle)
             ax.plot(-dx_, -dy_, "s", color=m["color"], markersize=2)
             ax.text(-dx_, -(dy_ + 10), "☋", color=m["color"],
            fontsize=8, ha="center", va="top", clip_on=True)

        #ax.plot(0, 0, "x", color="red", zorder=10)
        ax.add_patch(plt.Circle((0, 0), R_j, color="yellow", zorder=9))
        ax.text(0, -20, "Sun", color="yellow", ha="center")

    ax.set_xlim(-1.1*amax*(1+emax),1.1*amax*(1+emax))
    ax.set_ylim(-1.1*amax*(1+emax),1.1*amax*(1+emax))
    draw_ecliptic(ax, view_angle)
    ticks = np.linspace(-1.1*amax*(1+emax),1.1*amax*(1+emax),9)
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_yticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_xlabel("Distance in AU (1AU = 149,597,870 km)")
    ax.set_ylabel("Distance in AU (1AU = 149,597,870 km)")
    ax.grid(True, linestyle=':', linewidth=0.8, color='white', alpha=0.3)
    ax.set_axisbelow(True) 
    velocities = {}       
    for m in kb:   
        vx_v, vy_v, vz_v = body_velocity_3d(days, m)
        speed = np.sqrt(vx_v**2 + vy_v**2 + vz_v**2)
        velocities[m["name"]] = speed

    lines = ["Orbital Velocities (km/s)"]
    for name, v in sorted(velocities.items(), key=lambda x: -x[1]):
        lines.append(f"{name:8s} : {v:5.2f}")
        vel_text = "\n".join(lines)

    ax.text(0.98, 0.98,vel_text,transform=ax.transAxes,ha="right",va="top",fontsize=10,
            color="cyan", bbox=dict(facecolor="black",edgecolor="cyan",alpha=0.5,zorder=20))
    JD = 2451545.0 + days
    d, m, y = julian_to_gregorian(JD)
    if YEARS_TO_ANIMATE==0:
       ax.set_title(f"Kuiper Belt Objects at {d:02d}/{m:02d}/{y}\n"
                    f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n",
                    color="black")
       
       ax.text(
               0.02, 0.98,   
               "Aphelion = Q\n"
               "Perihelion = q\n"
               "Ascending Node = ☊\n"
               "Descending Node = ☋\n"
               "Above Ecliptic = Bright\n"
               "Below Ecliptic = Dim"
               ,
               transform=ax.transAxes,
               ha="left", va="top",
               fontsize=10, color="orange",
               bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                 )
    else:
        ax.set_title(f"Kuiper Belt Objects animation from {DD}/{MM}/{YY} to {d:02d}/{m:02d}/{y}\n"
                     f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n"
                     f"Simulation time = {YEARS_TO_ANIMATE:.0f} Years",
                     color="black")
        
        ax.text(
                0.02, 0.98,   
                "Aphelion = Q\n"
                "Perihelion = q\n"
                "Ascending Node = ☊\n"
                "Descending Node = ☋\n"
                "Above Ecliptic = Bright\n"
                "Below Ecliptic = Dim"
                ,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10, color="orange",
                bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                  )

def update_kb(frame):
    global cam_angle_deg, cam_angle, view_angle
    if DYNAMIC_VIEW == "y":
        frac = (frame / total_frames)**0.5
        frac = min(frac, 1.0)                 
        cam_angle_deg = 90 * (1 - frac)       
        cam_angle_deg = max(cam_angle_deg, 0) 
        cam_angle = np.deg2rad(cam_angle_deg)
        view_angle = np.pi/2 - cam_angle
    else:
         cam_angle = np.deg2rad(cam_angle_deg)
         view_angle = np.pi/2 - cam_angle
    
    draw_kb(ax_kb, start_days + frame*dt)
ani_kb = animation.FuncAnimation(fig_kb, update_kb,
                                        frames=total_frames, interval=17)
plt.show()


fig_kbo, ax_kbo = plt.subplots(figsize=(16,16))

def draw_kbo(ax, days):
    ax.clear()
    ax.set_facecolor("black")
    ax.set_aspect("equal")

    kbo = [m for m in moons if m["name"] in
             ["Eris","Sedna","2012 VP113","Leleākūhonua",
              "2023 KQ14","2004 VN112"]] 
    amax = max(m["a"] for m in kbo)
    emax = max(m["e"] for m in kbo)
    for m in kbo:
       X, Y, Z = compute_orbit_3d(days, m)
       xp, yp = project_2d(X, Y, Z, view_angle)
       ax.add_patch(plt.Circle((xp, yp), m["R"], color=m["color"]))
       ax.text(
           xp, yp,
           m["name"],
           color=m["color"],
           ha="center",
           va="center",
           fontsize=10,
           zorder=1,
           clip_on=True        
       )
       nu_ell = np.linspace(0, 2*pi, 400)
       Xo, Yo, Zo = [], [], []

       for nu in nu_ell:
           r = m["a"] * (1 - m["e"]**2) / (1 + m["e"] * np.cos(nu))
           x = r * np.cos(nu)
           y = r * np.sin(nu)

           ω = omega_at_time(m, days)
           cω, sω = np.cos(ω), np.sin(ω)
           x1 =  cω*x - sω*y
           y1 =  sω*x + cω*y

           ci, si = np.cos(m["i"]), np.sin(m["i"])
           x2 = x1
           y2 = ci*y1
           z2 = si*y1

           Omega_t = m["Omega"]
           cO, sO = np.cos(Omega_t), np.sin(Omega_t)
           Xo.append(cO*x2 - sO*y2)
           Yo.append(sO*x2 + cO*y2)
           Zo.append(z2)

       xx, yy = project_2d(-np.array(Xo), -np.array(Yo), np.array(Zo), view_angle)
       Z = np.array(Zo)
       xx_above, yy_above= xx.copy(), yy.copy()
       xx_below, yy_below = xx.copy(), yy.copy()
       xx_above[Z < 0], yy_above[Z < 0]  = np.nan, np.nan
       xx_below[Z >= 0], yy_below[Z >= 0] = np.nan, np.nan
       ls = LINESTYLE_MAP.get(m["type"])
       style = depth_style(Zo, ls["alpha"], ls["linewidth"])
       ax.plot(xx_above, yy_above,
       linestyle=ls["linestyle"],
       color=m["color"],
       alpha=style["above"]["alpha"],
       linewidth=style["above"]["linewidth"])
       ax.plot(xx_below, yy_below,
       linestyle=ls["linestyle"],
       color=m["color"],
       alpha=style["below"]["alpha"],
       linewidth=style["below"]["linewidth"])
       if m["e"] < 1:   
         (Xp, Yp, Zp), (Xa, Ya, Za) = peri_apo_3d(days, m)
         px, py = project_2d(Xp, Yp, Zp, view_angle)
         ax_, ay_ = project_2d(Xa, Ya, Za, view_angle)
         ax.plot(px, py, "o", color=m["color"], markersize=2, zorder=5)
         ax.text(px, py - 10*4, "q",
           color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
         ax.plot(ax_, ay_, "o", color=m["color"], markersize=2, zorder=5)
         ax.text(ax_, ay_ - 10*4, "Q",
           color=m["color"], ha="center", va="top", fontsize=9, clip_on=True)
         asc, desc = node_positions_3d(m, days)

         if asc is not None:
            ax_, ay_ = project_2d(*asc, view_angle)
            ax.plot(-ax_, -ay_, "s", color=m["color"], markersize=2)
            ax.text(-ax_, -(ay_ + 10*4), "☊", color=m["color"],
                    fontsize=8, ha="center", va="top", clip_on=True)

         if desc is not None:
            dx_, dy_ = project_2d(*desc, view_angle)
            ax.plot(-dx_, -dy_, "s", color=m["color"], markersize=2)
            ax.text(-dx_, -(dy_ + 10*4), "☋", color=m["color"],
           fontsize=8, ha="center", va="top", clip_on=True)

    #ax.plot(0, 0, "x", color="red", zorder=10)
    ax.add_patch(plt.Circle((0, 0), R_j, color="yellow", zorder=9))
    ax.text(0, -20, "Sun", color="yellow", ha="center")

    
    ax.set_xlim(-0.90*amax*(1+emax),0.90*amax*(1+emax))
    ax.set_ylim(-0.90*amax*(1+emax),0.90*amax*(1+emax))
    draw_ecliptic(ax, view_angle)
    ticks = np.linspace(-0.90*amax*(1+emax),0.90*amax*(1+emax),9)
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    ax.set_xticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_yticklabels([f"{t*c/149.59e6:.2f}" for t in ticks])
    ax.set_xlabel("Distance in AU (1AU = 149,597,870 km)")
    ax.set_ylabel("Distance in AU (1AU = 149,597,870 km)")
    ax.grid(True, linestyle=':', linewidth=0.8, color='white', alpha=0.3)
    ax.set_axisbelow(True) 
    velocities = {}       
    for m in kbo:   
        vx_v, vy_v, vz_v = body_velocity_3d(days, m)
        speed = np.sqrt(vx_v**2 + vy_v**2 + vz_v**2)
        velocities[m["name"]] = speed

    lines = ["Orbital Velocities (km/s)"]
    for name, v in sorted(velocities.items(), key=lambda x: -x[1]):
        lines.append(f"{name:8s} : {v:5.2f}")
        vel_text = "\n".join(lines)

    ax.text(0.98, 0.98,vel_text,transform=ax.transAxes,ha="right",va="top",fontsize=10,
            color="cyan", bbox=dict(facecolor="black",edgecolor="cyan",alpha=0.5,zorder=20))
    JD = 2451545.0 + days
    d, m, y = julian_to_gregorian(JD)
    if YEARS_TO_ANIMATE==0:
       ax.set_title(f"Extreme Trans-Neptunian Objects at {d:02d}/{m:02d}/{y}\n"
                    f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n",               
                    color="black")
       
       ax.text(
               0.02, 0.98,   
               "Aphelion = Q\n"
               "Perihelion = q\n"
               "Ascending Node = ☊\n"
               "Descending Node = ☋\n"
               "Above Ecliptic = Bright\n"
               "Below Ecliptic = Dim"
               ,
               transform=ax.transAxes,
               ha="left", va="top",
               fontsize=10, color="orange",
               bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                 )
    else:
        ax.set_title(f"Extreme Trans-Neptunian Objects animation from {DD}/{MM}/{YY} to {d:02d}/{m:02d}/{y}\n"
                     f"Viewing Angle = {cam_angle_deg:.0f} Degrees Above the Ecliptic\n"
                     f"Simulation time = {YEARS_TO_ANIMATE:.0f} Years",
                     color="black")
        
        ax.text(
                0.02, 0.98,   
                "Aphelion = Q\n"
                "Perihelion = q\n"
                "Ascending Node = ☊\n"
                "Descending Node = ☋\n"
                "Above Ecliptic = Bright\n"
                "Below Ecliptic = Dim"
                ,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=10, color="orange",
                bbox=dict(facecolor="black", edgecolor="orange", alpha=0.7)
                  )

def update_kbo(frame):
    global cam_angle_deg, cam_angle, view_angle
    if DYNAMIC_VIEW == "y":
        frac = (frame / total_frames)**0.5
        frac = min(frac, 1.0)                 
        cam_angle_deg = 90 * (1 - frac)       # 90° → 0°
        cam_angle_deg = max(cam_angle_deg, 0) # safety clamp
        cam_angle = np.deg2rad(cam_angle_deg)
        view_angle = np.pi/2 - cam_angle
    else:
         cam_angle = np.deg2rad(cam_angle_deg)
         view_angle = np.pi/2 - cam_angle
    
    draw_kbo(ax_kbo, start_days + frame*dt)
ani_kbo = animation.FuncAnimation(fig_kbo, update_kbo,
                                        frames=total_frames, interval=17)
plt.show()

# if FRAME_MODE == "h":
   
#    ani_zoom.save("C:/Users/aumva/Downloads/Heliocenter_Heliocentric.mp4", writer=writer)
# else:
#      FRAME_MODE == "b"
#      ani_zoom.save("C:/Users/aumva/Downloads/barycenter_barycentricnew.mp4", writer=writer)
     
#ani_inner.save("C:/Users/aumva/Downloads/inner_solarsystem1_halley.mp4", writer=writer)

# ani_outer.save("C:/Users/aumva/Downloads/outer_solarsystem1.mp4", writer=writer)

# ani_kb.save("C:/Users/aumva/Downloads/kuiper_belt1.mp4", writer=writer)

# ani_kbo.save("C:/Users/aumva/Downloads/outer_kuiper_belt1.mp4", writer=writer)



