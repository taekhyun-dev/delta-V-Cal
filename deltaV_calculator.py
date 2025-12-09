import math

# ============================
# 1. Physical constants
# ============================

G0 = 9.80665  # m/s^2

MU_SUN = 1.32712440018e20   # m^3/s^2
MU_EARTH = 3.986004418e14   # m^3/s^2
R_EARTH = 6371e3            # m
AU = 1.495978707e11         # m

R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# ============================
# 2. Hohmann Earth→Venus transfer
# ============================

def earth_venus_hohmann():
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
    ΔV from circular LEO to hyperbolic escape with given v_inf (m/s).
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv, v_circ, v_esc, v_peri_hyp


def venus_direct_entry_dv_required(leo_alt_km=300.0, margin_m_per_s=300.0):
    """
    KARI 시나리오: 금성 '대기 진입'만 고려 → VOI(궤도 포획) 없음.
    요구 ΔV = TVI (LEO->Venus transfer) + margin
    """
    hv = earth_venus_hohmann()
    dv_tvi, v_circ, v_esc, v_peri_hyp = dv_leo_to_vinf(hv["v_inf_earth"], leo_alt_km)

    dv_total = dv_tvi + margin_m_per_s
    return {
        "dv_tvi": dv_tvi,
        "dv_margin": margin_m_per_s,
        "dv_total": dv_total,
        "v_circ": v_circ,
        "v_esc": v_esc,
        "v_peri_hyp": v_peri_hyp,
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

def total_dv(cp_mass_kg):
    """
    Total ΔV available from Kick Stage + DSPS
    for a given Capsule Probe (CP) mass.

    cp_mass_kg : Capsule Probe mass [kg]

    Returns:
        total ΔV [m/s]
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

    # ---- KS burn ----
    # KS carries DSPS (struct + prop) + CP
    payload_for_ks = M_DSPS_STRUCT + M_DSPS_PROP + cp_mass_kg

    dv_ks = stage_delta_v(
        ISP_KS,
        M_KS_STRUCT,
        M_KS_PROP,
        payload_for_ks
    )

    # ---- DSPS burn ----
    # DSPS carries only CP
    dv_dsps = stage_delta_v(
        ISP_DSPS,
        M_DSPS_STRUCT,
        M_DSPS_PROP,
        cp_mass_kg
    )

    return dv_ks + dv_dsps


# ============================
# 4. Nuri + KS + DSPS + Capsule (KARI scenario)
# ============================

def nuri_venus_kari_scenario():
    # ---- Mission / environment ----
    LEO_ALT_KM = 300.0
    DV_MARGIN = 300.0  # m/s, mid-course & targeting margin

    # ---- Nuri LEO capacity ----
    M_NURI_LEO_LIMIT = 3200.0  # kg (3.2 t)

    # ---- Kick Stage (KS) from KARI Table 2 ----
    ISP_KS = 360.0
    M_KS_STRUCT = 450.0
    M_KS_PROP = 2050.0

    # ---- DSPS from KARI Table 2 ----
    ISP_DSPS = 290.0
    M_DSPS_STRUCT = 180.0
    M_DSPS_PROP = 120.0

    # ---- Capsule probe (design variable) ----
    M_CP = 400.0  # baseline from KARI [kg]

    # 1) ΔV requirement: TVI + margin (direct atmospheric entry)
    dv_req = venus_direct_entry_dv_required(
        leo_alt_km=LEO_ALT_KM,
        margin_m_per_s=DV_MARGIN
    )

    dv_tvi = dv_req["dv_tvi"]
    dv_total_req = dv_req["dv_total"]

    print("===== Venus direct entry ΔV budget (KARI-style) =====")
    print(f"LEO altitude         : {LEO_ALT_KM:.1f} km")
    print(f"ΔV (TVI, Earth dep.) : {dv_tvi/1000:.3f} km/s")
    print(f"ΔV margin            : {DV_MARGIN/1000:.3f} km/s")
    print(f"TOTAL ΔV required    : {dv_total_req/1000:.3f} km/s\n")

    # 2) Stage ΔV with given masses
    # KS carries: DSPS (struct+prop) + CP
    payload_for_KS = M_DSPS_STRUCT + M_DSPS_PROP + M_CP
    dv_ks = stage_delta_v(ISP_KS, M_KS_STRUCT, M_KS_PROP, payload_for_KS)

    # DSPS carries only the capsule
    payload_for_DSPS = M_CP
    dv_dsps = stage_delta_v(ISP_DSPS, M_DSPS_STRUCT, M_DSPS_PROP, payload_for_DSPS)

    dv_total_stages = dv_ks + dv_dsps

    print("===== Stage performance (from KARI Table 2) =====")
    print(f"KS ΔV       : {dv_ks/1000:.3f} km/s")
    print(f"DSPS ΔV     : {dv_dsps/1000:.3f} km/s")
    print(f"TOTAL ΔV    : {dv_total_stages/1000:.3f} km/s")

    # 3) Check Nuri LEO mass limit
    m_leo_total = (M_KS_STRUCT + M_KS_PROP +
                   M_DSPS_STRUCT + M_DSPS_PROP +
                   M_CP)

    print("\n===== LEO mass check =====")
    print(f"Total mass to LEO (KS+DSPS+CP): {m_leo_total:.1f} kg")
    print(f"Nuri LEO capacity             : {M_NURI_LEO_LIMIT:.1f} kg")

    # 4) Feasibility assessment
    print("\n===== Feasibility assessment =====")
    if m_leo_total > M_NURI_LEO_LIMIT:
        print("-> EXCEEDS Nuri LEO capacity: reduce stage/CP masses.")
    else:
        print("-> OK for Nuri LEO capacity.")

    if dv_total_stages < dv_total_req:
        print("-> ΔV INSUFFICIENT even for M_CP = 400 kg.")
    else:
        dv_margin = dv_total_stages - dv_total_req
        print("-> ΔV SUFFICIENT for direct Venus entry.")
        print(f"   ΔV margin: {dv_margin/1000:.3f} km/s")

    def find_cp_max():
        cp_max = 0.0
        for cp in range(0, 801):  # 0~800 kg까지 1 kg 단위로 살펴보기
            m_leo = 2500 + 300 + cp  # KS + DSPS + CP
            if m_leo > 3200:  # 누리 LEO 한계
                continue
            dv_total = total_dv(cp)  # KS+DSPS ΔV 합
            if dv_total >= dv_total_req:  # 요구 ΔV 만족?
                cp_max = cp
        print("Max CP mass (both constraints) =", cp_max, "kg")

    find_cp_max()



if __name__ == "__main__":
    nuri_venus_kari_scenario()
