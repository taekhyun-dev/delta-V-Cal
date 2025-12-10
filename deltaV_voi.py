import math

# ============================
# 1. Physical constants
# ============================

G0 = 9.80665  # m/s^2

MU_SUN = 1.32712440018e20   # m^3/s^2
MU_EARTH = 3.986004418e14   # m^3/s^2
R_EARTH = 6371e3            # m
AU = 1.495978707e11         # m

# Venus constants
R_VENUS = 6051.8e3          # m
MU_VENUS = 3.24859e14       # m^3/s^2

R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# ============================
# 2. Hohmann Earth→Venus transfer
# ============================

def earth_venus_hohmann():
    """Sun-centered Hohmann transfer between Earth (1 AU) and Venus (0.723 AU)."""
    r1 = R_EARTH_ORBIT
    r2 = R_VENUS_ORBIT
    a_trans = 0.5 * (r1 + r2)

    v_earth = math.sqrt(MU_SUN / r1)
    v_venus = math.sqrt(MU_SUN / r2)

    v_trans_earth = math.sqrt(MU_SUN * (2.0/r1 - 1.0/a_trans))
    v_trans_venus = math.sqrt(MU_SUN * (2.0/r2 - 1.0/a_trans))

    v_inf_earth = abs(v_trans_earth - v_earth)
    v_inf_venus = abs(v_trans_venus - v_venus)

    return {
        "v_earth": v_earth,
        "v_venus": v_venus,
        "v_trans_earth": v_trans_earth,
        "v_trans_venus": v_trans_venus,
        "v_inf_earth": v_inf_earth,
        "v_inf_venus": v_inf_venus,
    }


def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    """
    ΔV from circular LEO at leo_alt_km to a hyperbolic escape with v_inf (m/s).
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv, v_circ, v_esc, v_peri_hyp


def dv_venus_capture(v_inf_venus, venus_orbit_alt_km=300.0):
    """
    ΔV to capture into a circular Venus orbit at given altitude,
    coming in with hyperbolic excess speed v_inf_venus (m/s).
    """
    r_p = R_VENUS + venus_orbit_alt_km * 1e3

    v_circ = math.sqrt(MU_VENUS / r_p)
    v_esc = math.sqrt(2.0 * MU_VENUS / r_p)
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)

    dv_capture = v_peri_hyp - v_circ
    return dv_capture


def venus_orbit_dv_required(leo_alt_km=300.0,
                            venus_orbit_alt_km=300.0,
                            margin_m_per_s=300.0):
    """
    요구 ΔV = TVI (LEO -> Venus transfer) + VOI (Venus orbit insertion) + margin
    지구 복귀는 고려하지 않음.
    """
    hv = earth_venus_hohmann()

    # 1) TVI from LEO
    dv_tvi, _, _, _ = dv_leo_to_vinf(hv["v_inf_earth"], leo_alt_km)

    # 2) VOI at Venus
    dv_voi = dv_venus_capture(hv["v_inf_venus"], venus_orbit_alt_km)

    dv_total = dv_tvi + dv_voi + margin_m_per_s
    return {
        "dv_tvi": dv_tvi,
        "dv_voi": dv_voi,
        "dv_margin": margin_m_per_s,
        "dv_total": dv_total,
    }


# ============================
# 3. Rocket equation utilities
# ============================

def stage_delta_v(isp_s, m_struct, m_prop, m_payload):
    """
    Single stage ΔV with:
        isp_s    : specific impulse [s]
        m_struct : structural mass [kg]
        m_prop   : propellant mass [kg]
        m_payload: payload mass [kg] carried by this stage
    """
    if m_prop <= 0.0:
        return 0.0

    m0 = m_struct + m_prop + m_payload
    mf = m_struct + m_payload

    if m0 <= mf:
        return 0.0

    vex = isp_s * G0
    return vex * math.log(m0 / mf)


def total_dv_from_stages(cp_mass_kg):
    """
    Total ΔV available from Kick Stage + DSPS
    for a given Capsule Probe (CP) mass.
    """
    # ---- KARI Table 2 constants ----
    # Kick Stage (KS)
    ISP_KS = 360.0
    M_KS_STRUCT = 450.0
    M_KS_PROP = 2050.0

    # Deep Space Propulsion System (DSPS)
    ISP_DSPS = 290.0
    M_DSPS_STRUCT = 180.0
    M_DSPS_PROP = 120.0

    # KS carries DSPS (struct + prop) + CP
    payload_for_ks = M_DSPS_STRUCT + M_DSPS_PROP + cp_mass_kg
    dv_ks = stage_delta_v(
        ISP_KS, M_KS_STRUCT, M_KS_PROP, payload_for_ks
    )

    # DSPS carries only CP
    dv_dsps = stage_delta_v(
        ISP_DSPS, M_DSPS_STRUCT, M_DSPS_PROP, cp_mass_kg
    )

    return dv_ks + dv_dsps


# ============================
# 4. Main: Nuri + KS + DSPS + Capsule (orbit insertion only)
# ============================

def nuri_venus_orbit_scenario():
    # Mission settings
    LEO_ALT_KM = 300.0
    VENUS_ORBIT_ALT_KM = 300.0
    DV_MARGIN = 300.0  # m/s

    # Nuri LEO capacity
    M_NURI_LEO_LIMIT = 3200.0  # kg

    # KS / DSPS masses (from KARI)
    M_KS_STRUCT = 450.0
    M_KS_PROP = 2050.0
    M_DSPS_STRUCT = 180.0
    M_DSPS_PROP = 120.0

    # 1) Required ΔV for orbit insertion (no Earth return)
    dv_req = venus_orbit_dv_required(
        leo_alt_km=LEO_ALT_KM,
        venus_orbit_alt_km=VENUS_ORBIT_ALT_KM,
        margin_m_per_s=DV_MARGIN
    )
    dv_total_req = dv_req["dv_total"]

    print("===== Venus orbit insertion ΔV budget =====")
    print(f"LEO altitude                : {LEO_ALT_KM:.1f} km")
    print(f"Target Venus orbit altitude : {VENUS_ORBIT_ALT_KM:.1f} km")
    print(f"ΔV (TVI)                    : {dv_req['dv_tvi']/1000:.3f} km/s")
    print(f"ΔV (VOI)                    : {dv_req['dv_voi']/1000:.3f} km/s")
    print(f"ΔV margin                   : {dv_req['dv_margin']/1000:.3f} km/s")
    print(f"TOTAL ΔV required           : {dv_total_req/1000:.3f} km/s\n")

    # 2) Scan CP mass to find maximum feasible value
    cp_max = None
    for cp in range(0, 801):  # 0 ~ 800 kg
        # Mass to LEO
        m_leo_total = (M_KS_STRUCT + M_KS_PROP +
                       M_DSPS_STRUCT + M_DSPS_PROP +
                       cp)
        if m_leo_total > M_NURI_LEO_LIMIT:
            continue  # violates Nuri capacity

        dv_total_stages = total_dv_from_stages(cp)
        if dv_total_stages >= dv_total_req:
            cp_max = cp  # feasible; keep the largest

    # 3) Report feasibility
    print("===== Feasibility (with KS+DSPS from KARI) =====")
    if cp_max is None:
        print("-> No capsule mass satisfies both ΔV and Nuri capacity.")
        print("   (Even CP = 0 kg cannot reach Venus orbit with this stage design.)")
    else:
        print(f"-> Max feasible CP mass      : {cp_max:.1f} kg")
        dv_total_stages = total_dv_from_stages(cp_max)
        print(f"   ΔV from stages (CP={cp_max} kg): {dv_total_stages/1000:.3f} km/s")
        print(f"   ΔV margin                  : {(dv_total_stages - dv_total_req):.1f} m/s")


if __name__ == "__main__":
    nuri_venus_orbit_scenario()
