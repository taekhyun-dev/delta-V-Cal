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
R_VENUS_ORBIT = 0.723 * AU  # circular Venus orbit approx


# ============================
# 2. Hohmann Earth→Venus transfer
# ============================

def earth_venus_hohmann():
    r1 = R_EARTH_ORBIT
    r2 = R_VENUS_ORBIT
    a_trans = 0.5 * (r1 + r2)

    v_earth = math.sqrt(MU_SUN / r1)
    v_trans_earth = math.sqrt(MU_SUN * (2.0/r1 - 1.0/a_trans))

    v_inf_earth = abs(v_trans_earth - v_earth)
    return v_inf_earth


def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    """
    ΔV from circular LEO to hyperbolic escape with given v_inf (m/s).
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv  # [m/s]


def venus_entry_dv_required(leo_alt_km=300.0, margin_m_per_s=300.0):
    """
    금성 '대기 진입'만 고려하는 ΔV 요구:
    - LEO -> 금성 전이(TVI)
    - + margin (TCM 등)
    """
    v_inf_earth = earth_venus_hohmann()
    dv_tvi = dv_leo_to_vinf(v_inf_earth, leo_alt_km)
    dv_total = dv_tvi + margin_m_per_s
    return dv_total, dv_tvi


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

    cp_mass_kg : Capsule Probe mass [kg]
    Returns:
        total ΔV [m/s]
    """
    # Kick Stage (KS) from KARI
    ISP_KS = 360.0
    M_KS_STRUCT = 450.0
    M_KS_PROP = 2050.0

    # Deep Space Propulsion System (DSPS)
    ISP_DSPS = 290.0
    M_DSPS_STRUCT = 180.0
    M_DSPS_PROP = 120.0

    # KS burn: carries DSPS (struct+prop) + CP
    payload_for_ks = M_DSPS_STRUCT + M_DSPS_PROP + cp_mass_kg
    dv_ks = stage_delta_v(ISP_KS, M_KS_STRUCT, M_KS_PROP, payload_for_ks)

    # DSPS burn: carries only CP
    payload_for_dsps = cp_mass_kg
    dv_dsps = stage_delta_v(ISP_DSPS, M_DSPS_STRUCT, M_DSPS_PROP, payload_for_dsps)

    return dv_ks + dv_dsps


# ============================
# 4. Solve for max CP mass (entry case)
# ============================

def find_cp_max_for_entry():
    # Nuri LEO capacity (KARI 가정)
    M_NURI_LEO_LIMIT = 3200.0  # kg

    # KS + DSPS fixed masses
    M_KS_TOTAL = 2500.0  # 2050 + 450
    M_DSPS_TOTAL = 300.0 # 120 + 180

    # 1) ΔV 요구 계산
    dv_req_total, dv_tvi = venus_entry_dv_required(
        leo_alt_km=300.0,
        margin_m_per_s=300.0
    )

    print("===== Venus atmospheric entry case =====")
    print(f"Required ΔV (TVI + margin) : {dv_req_total/1000:.3f} km/s")
    print(f"  - TVI only                : {dv_tvi/1000:.3f} km/s\n")

    # 2) CP 스윕하면서 두 조건을 만족하는 최대값 찾기
    cp_max = 0.0
    cp_max_dv_only = 0.0

    for cp in range(0, 1001):  # 0 ~ 1000 kg까지 1kg 단위 탐색
        cp_mass = float(cp)

        # Nuri LEO 질량 제한 체크
        m_leo = M_KS_TOTAL + M_DSPS_TOTAL + cp_mass
        if m_leo > M_NURI_LEO_LIMIT:
            continue  # 발사체 한계 초과

        dv_total = total_dv_from_stages(cp_mass)

        # ΔV 요구 만족 여부
        if dv_total >= dv_req_total:
            cp_max = cp_mass  # LEO + ΔV 둘 다 만족하는 최대 CP

        # “ΔV 조건만” 만족하는 최대 CP도 따로 기록
        if dv_total >= dv_req_total:
            cp_max_dv_only = cp_mass

    print(f"Max CP mass (ΔV + Nuri LEO limit) : {cp_max:.1f} kg")
    print(f"(for reference) Max CP mass (ΔV only, ignoring launcher) : {cp_max_dv_only:.1f} kg")

    # 참고: CP=0, CP=cp_max에서 실제 ΔV 출력
    dv_cp0 = total_dv_from_stages(0.0) / 1000.0
    dv_cpmax = total_dv_from_stages(cp_max) / 1000.0
    print(f"Total stage ΔV with CP=0 kg       : {dv_cp0:.3f} km/s")
    print(f"Total stage ΔV with CP={cp_max:.1f} kg : {dv_cpmax:.3f} km/s")


if __name__ == "__main__":
    find_cp_max_for_entry()
