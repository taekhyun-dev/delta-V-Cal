import math
import csv

# ============================
# 1. Physical constants
# ============================

G0 = 9.80665  # m/s^2

# Gravitational parameters [m^3/s^2]
MU_SUN   = 1.32712440018e20
MU_EARTH = 3.986004418e14
MU_VENUS = 3.24859e14

# Planet radii [m]
R_EARTH = 6371e3
R_VENUS = 6051.8e3

# Orbits around the Sun
AU = 1.495978707e11
R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# ============================
# 2. Heliocentric Earth↔Venus transfer
# ============================

def earth_venus_hohmann():
    """
    Compute a simple heliocentric Hohmann transfer from Earth (1 AU)
    to Venus (~0.723 AU).

    Returns:
        dict with v_inf at Earth and Venus, etc.
    """
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


# ============================
# 3. ΔV from LEO to v_inf (for KS)
# ============================

def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    """
    ΔV from circular LEO at altitude leo_alt_km to a hyperbolic escape
    with excess velocity v_inf (m/s) relative to Earth.

    This is what the Kick Stage (KS) must provide for TVI.
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc  = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv, v_circ, v_esc, v_peri_hyp


# ============================
# 4. ΔV for Venus Orbit Insertion (for DSPS)
# ============================

def dv_venus_orbit_insertion(v_inf_venus, venus_orbit_alt_km=300.0):
    """
    ΔV to capture from a hyperbolic approach with excess speed v_inf_venus
    into a circular orbit around Venus at altitude venus_orbit_alt_km.

    This is what the DSPS must provide for VOI (not atmospheric entry).
    """
    r_p = R_VENUS + venus_orbit_alt_km * 1e3

    v_circ = math.sqrt(MU_VENUS / r_p)          # circular speed at target orbit
    v_esc  = math.sqrt(2.0 * MU_VENUS / r_p)    # local escape speed
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)

    dv_capture = v_peri_hyp - v_circ
    return dv_capture, v_circ, v_peri_hyp


# ============================
# 5. Rocket equation
# ============================

def stage_delta_v(isp_s, m_struct, m_prop, m_payload):
    """
    Single stage ΔV (m/s) for given:
      isp_s    : specific impulse [s]
      m_struct : structural mass [kg]
      m_prop   : propellant mass [kg]
      m_payload: payload mass carried by this stage [kg]
    """
    if m_prop <= 0.0:
        return 0.0

    m0 = m_struct + m_prop + m_payload
    mf = m_struct + m_payload

    if m0 <= mf:
        return 0.0

    vex = isp_s * G0
    return vex * math.log(m0 / mf)


# ============================
# 6. ΔV requirements for KS and DSPS in our scenario
# ============================

def dv_requirements_for_orbit_insertion(
    leo_alt_km=300.0,
    venus_orbit_alt_km=300.0,
    dsps_margin_m_per_s=300.0
):
    """
    Our scenario:
      - KS: LEO -> Venus transfer (TVI)
      - DSPS: TCM + VOI into Venus orbit

    dsps_margin_m_per_s: extra ΔV margin (for TCM, dispersions, etc.)
                          added on top of VOI requirement.
    """
    hv = earth_venus_hohmann()

    # KS requirement: TVI from LEO
    dv_tvi, v_circ_leo, v_esc_leo, v_peri_hyp_leo = dv_leo_to_vinf(
        hv["v_inf_earth"], leo_alt_km
    )

    # DSPS requirement: VOI (Venus orbit insertion) + margin
    dv_voi, v_circ_venus, v_peri_venus = dv_venus_orbit_insertion(
        hv["v_inf_venus"], venus_orbit_alt_km
    )
    dv_dsps_req = dv_voi + dsps_margin_m_per_s

    return {
        "dv_ks_req": dv_tvi,
        "dv_dsps_req": dv_dsps_req,
        "dv_voi": dv_voi,
        "v_inf_earth": hv["v_inf_earth"],
        "v_inf_venus": hv["v_inf_venus"],
        "v_circ_leo": v_circ_leo,
        "v_esc_leo": v_esc_leo,
        "v_peri_hyp_leo": v_peri_hyp_leo,
        "v_circ_venus": v_circ_venus,
        "v_peri_venus": v_peri_venus,
    }


# ============================
# 7. CP 질량 스윕해서 최대 CP 찾기
# ============================

def find_cp_max_for_orbit_mission(
    # Nuri LEO capability
    nuri_leo_limit_kg=3200.0,
    leo_alt_km=300.0,
    venus_orbit_alt_km=300.0,
    dsps_margin_m_per_s=300.0,
    # KS (Kick Stage) design (default: KARI)
    isp_ks_s=360.0,
    m_ks_struct_kg=450.0,
    m_ks_prop_kg=2050.0,
    # DSPS design (default: KARI)
    isp_dsps_s=290.0,
    m_dsps_struct_kg=180.0,
    m_dsps_prop_kg=120.0,
    cp_search_step_kg=1.0
):
    """
    스펙을 넣으면:
      - KS/DSPS가 각자 담당해야 할 ΔV를 계산하고
      - CP 질량을 0 ~ LEO 한계까지 스윕해서
      - ΔV 조건 + LEO 질량 조건을 모두 만족하는 CP_max를 찾는다.
    """

    # 1) ΔV 요구량 계산 (우리 시나리오용)
    dv_req = dv_requirements_for_orbit_insertion(
        leo_alt_km=leo_alt_km,
        venus_orbit_alt_km=venus_orbit_alt_km,
        dsps_margin_m_per_s=dsps_margin_m_per_s
    )

    dv_ks_req   = dv_req["dv_ks_req"]
    dv_dsps_req = dv_req["dv_dsps_req"]

    # 2) LEO 질량 제약에서 나오는 CP 상한
    base_mass = (m_ks_struct_kg + m_ks_prop_kg +
                 m_dsps_struct_kg + m_dsps_prop_kg)
    cp_max_by_leo = max(0.0, nuri_leo_limit_kg - base_mass)
    print(f"Debug: base_mass = {base_mass:.1f} kg, cp_max_by_leo = {cp_max_by_leo:.1f} kg".format())

    cp_max_feasible = -1.0  # -1이면 "불가능" 표시

    simulation_logs = []

    # 3) CP 질량 스윕
    cp = 0.0
    while cp <= cp_max_by_leo + 1e-6:
        # KS 단계에서 payload = DSPS(구조+추진제) + CP
        payload_for_ks = m_dsps_struct_kg + m_dsps_prop_kg + cp
        dv_ks_avail = stage_delta_v(
            isp_ks_s,
            m_ks_struct_kg,
            m_ks_prop_kg,
            payload_for_ks
        )

        # DSPS 단계에서 payload = CP
        payload_for_dsps = cp
        dv_dsps_avail = stage_delta_v(
            isp_dsps_s,
            m_dsps_struct_kg,
            m_dsps_prop_kg,
            payload_for_dsps
        )
        ks_pass = dv_ks_avail >= dv_ks_req
        dsps_pass = dv_dsps_avail >= dv_dsps_req
        total_success = ks_pass and dsps_pass
        # ΔV 조건: 각자 담당 ΔV를 만족해야 함
        if dv_ks_avail >= dv_ks_req and dv_dsps_avail >= dv_dsps_req:
            cp_max_feasible = cp  # 이 CP까지는 충분
        # else: 이 CP에서는 실패 → 더 큰 CP는 보통 더 나쁘지만,
        #       우리는 끝까지 스윕해서 "최대"를 찾는다.

        log_entry = {
            "cp_kg": cp,
            "ks_avail": dv_ks_avail,
            "ks_req": dv_ks_req,
            "ks_pass": ks_pass,
            "dsps_avail": dv_dsps_avail,
            "dsps_req": dv_dsps_req,
            "dsps_pass": dsps_pass,
            "total_success": total_success
        }
        simulation_logs.append(log_entry)
        cp += cp_search_step_kg

        # 3. CSV 파일로 상세 로그 저장
    file_name = "mission_simulation_log.csv"
    
    # 'utf-8-sig'는 엑셀에서 한글이나 특수문자(Δ)가 깨지지 않게 해줍니다.
    with open(file_name, mode='w', newline='', encoding='utf-8-sig') as file:
        writer = csv.writer(file)
        
        # (1) 헤더(제목 행) 작성
        header = [
            "CP_Weight(kg)", 
            "KS_Avail(m/s)", "KS_Req(m/s)", "KS_Pass", 
            "DSPS_Avail(m/s)", "DSPS_Req(m/s)", "DSPS_Pass", 
            "Mission_Success"
        ]
        writer.writerow(header)

        # (2) 데이터 한 줄씩 작성
        for log in simulation_logs:
            row = [
                log['cp_kg'],
                round(log['ks_avail'], 2),  # 소수점 2자리 반올림
                round(log['ks_req'], 2),
                "PASS" if log['ks_pass'] else "FAIL",
                round(log['dsps_avail'], 2),
                round(log['dsps_req'], 2),
                "PASS" if log['dsps_pass'] else "FAIL",
                "SUCCESS" if log['total_success'] else "FAIL"
            ]
            writer.writerow(row)

    return {
        "cp_max_kg": cp_max_feasible,
        "dv_ks_req": dv_ks_req,
        "dv_dsps_req": dv_dsps_req,
        "cp_max_by_leo": cp_max_by_leo,
    }


# ============================
# 8. Example main
# ============================

def main():
    # 기본값: KARI 설계, Nuri LEO limit = 3200 kg, margin = 300 m/s
    result = find_cp_max_for_orbit_mission(
        nuri_leo_limit_kg=3200.0,
        leo_alt_km=300.0,
        venus_orbit_alt_km=300.0,
        dsps_margin_m_per_s=300.0,
        isp_ks_s=360.0,
        m_ks_struct_kg=450.0,
        m_ks_prop_kg=2050.0,
        isp_dsps_s=290.0,
        m_dsps_struct_kg=180.0,
        m_dsps_prop_kg=120.0,
        cp_search_step_kg=1.0
    )

    print("===== Venus orbit mission with Nuri + KS + DSPS =====")
    print(f"Required ΔV (KS, TVI)       : {result['dv_ks_req']/1000:.3f} km/s")
    print(f"Required ΔV (DSPS, VOI+TCM) : {result['dv_dsps_req']/1000:.3f} km/s")
    print(f"CP max by LEO mass constraint : {result['cp_max_by_leo']:.1f} kg")
    print(f"CP max satisfying ΔV + LEO   : {result['cp_max_kg']:.1f} kg")

    if result["cp_max_kg"] < 0:
        print("=> Even CP = 0 kg cannot satisfy ΔV constraints with this design.")
    else:
        print("=> Up to this CP, KS and DSPS both have enough ΔV,")
        print("   and total mass stays within Nuri's LEO capability.")


if __name__ == "__main__":
    main()
