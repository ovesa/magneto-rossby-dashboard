# Visualize the sectoral dispersion curves for HD Rossby 
# & fast/slow magneto-Rossby (Eq. 44) + Alfvén, as a 
# function of azimuthal order m and magnetic field B.
# This is a Streamlit app that allows interactive exploration 
# of the parameter space.

###########################################
###########################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
import matplotlib

###########################################
################ Functions ################
###########################################

# magnetic permeability of free space in SI units (N A^-2)
# set the coupling strenght between magnetic fields and fluid motions.
mu0 = 4.0 * np.pi * 1e-7

def vA_from_B_gauss(B_G, rho, mu0=mu0):
    """Compute the Alfven speed from a magnetic field strength
    specified in Gauss. The Alfven speed is defined as
        
        vA = B / sqrt(mu0 * rho), 

    where B is the magnetic field in Tesla, mu0 is the permeability 
    of free space, and rho is the mass density.

    This function assumes the input magnetic field is provided in 
    Gauss and internally converts it to Tesla (1 G = 1e-4 T).

    Args:
        B_G (float or numpy.ndarray): Magnetic field strength in Gauss
        rho (float): Mass density in kg/m^3
        mu0 (float): Permeability of free space in SI units (N A^-2). 

    Returns:
        alfven_speed (float or numpy.ndarray): Alfven speed in m/s. 
                                The return type matches the type of B_G.
    """    
    B_T = B_G * 1e-4
    alfven_speed = B_T / np.sqrt(mu0 * rho)
    return alfven_speed

def nu_alfven_nHz(m, vA, R0):
    """Compute the cycle Alfven frequency for the sectoral modes.
    This function evalues the toroidal Alfven frequency given by 
    the dispersion relation
    
        omega = ± (m/R0)*vA,
        
    where m is the azimuthal order, R0 is the radius of the sphere,
    and vA is the Alfven speed. The frequency is then converted to
    cyclic frequency in nHz by dividing by 2π and multiplying by 1e9:
    
        nu_nHz = omega / (2π) * 1e9.
        
    Args:
        m (float or numpy.ndarray):Azimuthal wavenumber (dimensionless).
        vA (float): Alfvén speed in m s^-1.
        R0 (float): Characteristic radius in meters (e.g., solar radius).

    Returns:
        alfven_branches_nHz (float or numpy.ndarray): Cyclic Alfven frequency 
            in nHz. The return type matches the type of m.
            
    Notes:
        - The returned frequency corresponds to the magnitude of the 
            Alfven frequency; the dispersion relation has both positive 
            and negative branches, which are symmetric about zero. The sign
            must be applied externally if both prograde and retograde branches 
            are desired.
        - This function assumes a sectoral mode structure where the
            toroidal wavenumber is equal.
    """    
    alfven_branches_nHz = (m * vA / R0) / (2.0 * np.pi) * 1e9
    return alfven_branches_nHz


def solve_dispersion_eqn(s, n, alpha):
    """Solve the quadratic magneto-Rossby wave spherical dispersion relation 
    given by Eq. (44) in Zaqarashvili, Oliver, Ballester, and Shergelashvili (2007).
    This function sovles the quadratic equation for x = lambda/s:
    
        n(n+1) x^2 + x + alpha^2 [2 - n(n+1)] = 0,
        
    where lambda is the dimensionless eigenvalue,n is the poloidal wavenumber, 
    sis the toroidal wavenumber, and alpha is the dimensionless magnetic 
    parameter:
    
        alpha = vA / (2 * Omega0 * R0).
        
    The two returned roots correspond to the fast and slow magneto-Rossby 
    branches. Branch classification (fast vs slow) is handled externally.
    
    Args:
        s (float or numpy.ndarray): Toroidal (azimuthal) mode index. 
                                    For sectoral modes, s = m.
        n (float or numpy.ndarray): Poloidal mode index. 
                                    For sectoral modes, n = s = m.
        alpha (float): Dimensionless magnetic parameter defined as
                        vA / (2 Omega R0).

    Returns:
        lam1, lam2 tuple[numpy.ndarray, numpy.ndarray]: The two roots of the dispersion 
                        equation. The return type matches the type of s and n.
                        
    Notes:
        - In the limit alpha -> 0, one root approaches the hydrodynamic 
            Rossby solution: lambda/s = -1/[n(n+1)].
    """
    A = n * (n + 1.0)
    B = 1.0
    C = alpha**2 * (2.0 - n * (n + 1.0))

    disc = B * B - 4.0 * A * C
    # In typical regimes, disc >= 0. Guard small negative due to floating error:
    disc = np.maximum(disc, 0.0)
    sqrt_disc = np.sqrt(disc)

    x1 = (-B + sqrt_disc) / (2.0 * A)
    x2 = (-B - sqrt_disc) / (2.0 * A)

    lam1 = s * x1
    lam2 = s * x2
    return lam1, lam2

def lambda_to_nu_nHz(lam, nu_Omega_nHz):
    """Convert the dimensionless eigenvalue lambda to cyclic frequency in nHz.

    In the spherical magneto-Rossby formulation, the mode frequency is
    expressed in dimensionless form as

        omega / (2 Omega) = lambda,

    which implies

        nu_mode = 2 * nu_Omega * lambda,

    where nu_Omega is the cyclic rotation frequency of the background
    sphere and nu_mode is the cyclic mode frequency.

    Args:
        lam (float or numpy.ndarray): Dimensionless eigenvalue lambda 
                    obtained from the dispersion relation.
        nu_Omega_nHz (float): Background rotation frequency in nHz.

    Returns:
        lam_nHz (float or numpy.ndarray): Mode cyclic frequency in nHz. 
                The return type matches the type of lam.
                
    Notes:
        - The sign of lambda determines the prograde (positive) vs 
            retrograde (negative) nature of the mode in the rotating frame.

    """    
    lam_nHz = 2.0 * nu_Omega_nHz * lam
    return lam_nHz

def hd_rossby_nu_nHz_sectoral(m, nu_Omega_nHz):
    """Compute the hydrodynamic sectoral Rossby mode frequency in nHz.
    This function evaluates the alpha = 0 limit of the spherical
    magneto-Rossby dispersion relation (Eq. 44), corresponding to
    the purely hydrodynamic case. In this limit,

        lambda / s = -1 / [n(n+1)],

    and for sectoral modes (n = s = m),

        lambda = -m / [m(m+1)].

    The cyclic frequency is then obtained using

        nu = 2 * nu_Omega * lambda.

    Args:
        m (float or numpy.ndarray): Azimuthal (sectoral) wavenumber, 
                                    where n = s = m.
        nu_Omega_nHz (float): Background rotation frequency in nHz.

    Returns:
        res (float or numpy.ndarray): Cyclic Rossby mode frequency in nHz. 
                The return type matches the type of m.
                
    Notes:
        - The resulting frequency is negative, corresponding to retrograde 
            propagation in the rotating frame. 
    """    
    A = m * (m + 1.0)
    x_hd = -1.0 / A
    lam_hd = m * x_hd
    res = lambda_to_nu_nHz(lam_hd, nu_Omega_nHz)
    return res

###########################################
############### Streamlit UI ##############
###########################################

st.set_page_config(page_title="Magneto-Rossby Dispersion Equation", layout="wide")
st.markdown("""
<style>

/* Sidebar section headers */
section[data-testid="stSidebar"] h2 {
    font-size: 24px !important;
}

/* Sidebar number_input / selectbox / text_input labels */
section[data-testid="stSidebar"] label p {
    font-size: 18px !important;
}

/* Sidebar divider spacing */
section[data-testid="stSidebar"] hr {
    margin-top: 1.2rem;
    margin-bottom: 1.2rem;
}

</style>
""", unsafe_allow_html=True)

st.title("Dispersion Equation for Magneto-Rossby and Alfvén Waves")

with st.sidebar:
    st.header("Background parameters")

    nu_Omega_nHz = st.number_input(r"Rotation $\Omega$ [nHz]", value=456.0, step=1.0, format="%.1f")
    R0 = st.number_input(r"Solar radius $R_{\odot}$ [m]", value=6.96e8, step=1.0e7, format="%.3e")
    rho = st.number_input(r"Density $\rho$ [kg m$^{-3}$]", value=1e-2, step=1e-2, format="%.3e")

    st.divider()
    st.header("Mode range")
    m_min = st.number_input(r"azimuthal order $m$ min", min_value=1, max_value=500, value=1, step=1)
    m_max = st.number_input(r"azimuthal orde $m$ max", min_value=1, max_value=500, value=30, step=1)
    if m_max < m_min:
        st.warning("m max < m min; swapping.")
        m_min, m_max = m_max, m_min

    st.divider()
    st.header("Magnetic field(s) [G]")

    preset = st.selectbox("B presets", ["Custom", "0,10,30,50", "0,5,10,20,50,100", "0..100 step 10"])
    if preset == "0,10,30,50":
        B_list_G = [0, 10, 30, 50]
    elif preset == "0,5,10,20,50,100":
        B_list_G = [0, 5, 10, 20, 50, 100]
    elif preset == "0..100 step 10":
        B_list_G = list(range(0, 101, 10))
    else:
        B_text = st.text_input("Custom B list in Gauss (comma-separated)", value="0, 10, 30, 50")
        try:
            B_list_G = [float(x.strip()) for x in B_text.split(",") if x.strip() != ""]
            if len(B_list_G) == 0:
                B_list_G = [0.0]
        except Exception:
            st.error("Could not parse B list. Using [0, 10, 30, 50].")
            B_list_G = [0.0, 10.0, 30.0, 50.0]

    st.divider()
    st.markdown(
    "<h2 style='font-size:28px;'>Plot toggles</h2>",
    unsafe_allow_html=True
)
    show_hd = st.checkbox("Show HD Rossby baseline", value=True)
    show_fast_slow = st.checkbox("Show fast/slow magneto-Rossby branches", value=True)
    show_alfven = st.checkbox("Show Alfvén branch", value=True)

    st.divider()
    st.header("Axes")
    xlim0, xlim1 = st.slider("$m$ limit", min_value=0, max_value=max(60, int(m_max)), value=(0, int(m_max)))
    ylim0, ylim1 = st.slider(r"$\nu$ limit (nHz)", min_value=-5000, max_value=5000, value=(-600, 300))

###########################################
############### Compute ###################
###########################################

nu_Omega_Hz = nu_Omega_nHz * 1e-9
Omega0 = 2.0 * np.pi * nu_Omega_Hz  # rad/s

# Azimuthal order array for sectoral modes (n = s = m)
m = np.arange(m_min, m_max + 1, dtype=float)
s = m.copy()
n = m.copy()

# HD baseline + HD x for classification
x_hd = -1.0 / (m * (m + 1.0))
nu_hd = hd_rossby_nu_nHz_sectoral(m, nu_Omega_nHz)

# Build a tidy table of all curves (useful for export)
rows = []
if show_hd:
    for mi, nui in zip(m, nu_hd):
        rows.append({"curve": "HD Rossby", "B_G": 0.0, "m": mi, "nu_nHz": nui})

for B_G in B_list_G:
    if B_G == 0:
        continue

    vA = vA_from_B_gauss(B_G, rho)
    alpha = vA / (2.0 * Omega0 * R0)

    lam1, lam2 = solve_dispersion_eqn(s, n, alpha)

    x1 = lam1 / s
    x2 = lam2 / s
    fast_mask = np.abs(x1 - x_hd) <= np.abs(x2 - x_hd)

    lam_fast = np.where(fast_mask, lam1, lam2)
    lam_slow = np.where(fast_mask, lam2, lam1)

    nu_fast = lambda_to_nu_nHz(lam_fast, nu_Omega_nHz)
    nu_slow = lambda_to_nu_nHz(lam_slow, nu_Omega_nHz)

    if show_fast_slow:
        for mi, nui in zip(m, nu_fast):
            rows.append({"curve": "Fast magneto-Rossby", "B_G": float(B_G), "m": mi, "nu_nHz": nui})
        for mi, nui in zip(m, nu_slow):
            rows.append({"curve": "Slow magneto-Rossby", "B_G": float(B_G), "m": mi, "nu_nHz": nui})

    if show_alfven:
        nuA = nu_alfven_nHz(m, vA, R0)
        for mi, nui in zip(m, +nuA):
            rows.append({"curve": "Alfvén +", "B_G": float(B_G), "m": mi, "nu_nHz": nui})
        for mi, nui in zip(m, -nuA):
            rows.append({"curve": "Alfvén −", "B_G": float(B_G), "m": mi, "nu_nHz": nui})

df = pd.DataFrame(rows)

###########################################
######### Plotting with Streamlit #########
###########################################

fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(111)
ax.axhline(0, linewidth=1, linestyle="-.", color="k")

B_nonzero = sorted(set(float(b) for b in B_list_G if float(b) != 0.0))

cmap = matplotlib.colormaps.get_cmap("Set1")
max_colors = 15
idx = np.linspace(0, max_colors - 1, num=len(B_nonzero), dtype=int)
color_list = [cmap(i) for i in idx]

B_color = {B: color_list[i] for i, B in enumerate(B_nonzero)}
B_color[0.0] = "black"

def color_for_B(B_G):
    return B_color.get(float(B_G), "gray")

if show_hd:
    ax.plot(m, nu_hd, linewidth=2.5, color="black",
            label=r"HD Rossby")

for B_G in B_list_G:
    if B_G == 0:
        continue

    c = color_for_B(B_G)

    vA = vA_from_B_gauss(B_G, rho)
    alpha = vA / (2.0 * Omega0 * R0)

    lam1, lam2 = solve_dispersion_eqn(s, n, alpha)

    x1 = lam1 / s
    x2 = lam2 / s
    fast_mask = np.abs(x1 - x_hd) <= np.abs(x2 - x_hd)
    lam_fast = np.where(fast_mask, lam1, lam2)
    lam_slow = np.where(fast_mask, lam2, lam1)

    nu_fast = lambda_to_nu_nHz(lam_fast, nu_Omega_nHz)
    nu_slow = lambda_to_nu_nHz(lam_slow, nu_Omega_nHz)

    if show_fast_slow:
        ax.plot(m, nu_fast, linestyle="dotted", linewidth=3, color=c,
                label=f"Fast [B={B_G:g} G]")
        ax.plot(m, nu_slow, linestyle="-", linewidth=3, color=c,
                label=f"Slow [B={B_G:g} G]")

    if show_alfven:
        nuA = nu_alfven_nHz(m, vA, R0)
        ax.plot(m, +nuA, linestyle="dashdot", linewidth=2, color=c,
                label=f"Alfvén [B={B_G:g} G]")
        ax.plot(m, -nuA, linestyle="dashdot", linewidth=2, color=c,
                label="_nolegend_")  # avoids duplicate legend entry

ax.set_xlabel(r"$m$", fontsize=18)
ax.set_ylabel(r"$\nu$ [nHz]", fontsize=18)
ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.2)
ax.tick_params(axis="both", which="minor", labelsize=16, length=3, width=1.0)
ax.minorticks_on()

ax.set_xlim(xlim0, xlim1)
ax.set_ylim(ylim0, ylim1)
ax.grid(True, alpha=0.3)

ax.set_title(
    "Dispersion Equation for Sectoral Magneto-Rossby and Alfvén Waves\n"
    + r"$\rho$ = " + f"{rho:.3e}" + r" kg m$^{{-3}}$, $\Omega/(2\pi)$ = "
    + f"{nu_Omega_nHz:.1f}" + r" nHz",
    fontsize=20,
)
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
          fontsize=14, frameon=False)

fig.subplots_adjust(right=0.75)

st.pyplot(fig, clear_figure=True)

st.subheader("Data")
st.dataframe(df, use_container_width=True, height=420)

csv = df.to_csv(index=False).encode("utf-8")
st.download_button(
    "Download CSV",
    data=csv,
    file_name="magnetic_dispersion_curves.csv",
    mime="text/csv",
)
