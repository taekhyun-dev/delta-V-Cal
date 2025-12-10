import math

# =========================================
# 1. Physical constants
# =========================================

G0 = 9.80665  # m/s^2

MU_SUN = 1.32712440018e20   # m^3/s^2
MU_EARTH = 3.986004418e14   # m^3/s^2
MU_VENUS = 3.24859e14       # m^3/s^2

R_EARTH = 6371e3            # m
R_VENUS = 6051.8e3          # m
AU = 1.495978707e11         # m

R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# =========================================
# 2. Hohmann Earth→Venus transfer utilities
# =========================================

def earth_venus_hohmann():
    """Heliocentric Hohmann transfer Earth↔Venus (circular orbits)."""
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


def dv_venus_capture(v_inf_venus, venus_orbit_alt_km=300.0):
    """
    ΔV for Venus orbit insertion from a hyperbolic approach
    with given v_inf_venus (m/s) to circular orbit at venus_orbit_alt_km.
    """
    r_p = R_VENUS + venus_orbit_alt_km * 1e3

    v_circ = math.sqrt(MU_VENUS / r_p)
    v_esc = math.sqrt(2.0 * MU_VENUS / r_p)
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)

    dv_capture = v_peri_hyp - v_circ
    return dv_capture


# =========================================
# 3. ΔV requirements (direct_entry vs orbit_insertion)
# =========================================

def dv_required(
    mode: str = "direct_entry",
    leo_alt_km: float = 300.0,
    venus_orbit_alt_km: float = 300.0,
    margin_m_per_s: float = 300.0,
):
    """
    mode:
        - "direct_entry"   : TVI + margin (KARI 스타일 대기진입)
        - "orbit_insertion": TVI + VOI + margin (금성 궤도 포획)

    Returns:
        dv_total_req [m/s]
    """
    hv = earth_venus_hohmann()

    # TVI (LEO -> Venus transfer)
    dv_tvi, _, _, _ = dv_leo_to_vinf(hv["v_inf_earth"], leo_alt_km)

    if mode == "direct_entry":
        dv_voi = 0.0
    elif mode == "orbit_insertion":
        dv_voi = dv_venus_capture(hv["v_inf_venus"], venus_orbit_alt_km)
    else:
        raise ValueError("mode must be 'direct_entry' or 'orbit_insertion'")

    dv_total = dv_tvi + dv_voi + margin_m_per_s

    return {
        "dv_total": dv_total,
        "dv_tvi": dv_tvi,
        "dv_voi": dv_voi,
        "dv_margin": margin_m_per_s,
    }


# =========================================
# 4. Rocket equation utilities
# =========================================

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


def total_dv_for_cp(
    cp_mass_kg,
    ks_isp_s,
    ks_struct_kg,
    ks_prop_kg,
    dsps_isp_s,
    dsps_struct_kg,
    dsps_prop_kg,
):
    """
    Total ΔV available from Kick Stage + DSPS for a given CP mass.

    cp_mass_kg : Capsule Probe mass [kg]
    """

    # KS carries: DSPS (struct + prop) + CP
    payload_for_ks = dsps_struct_kg + dsps_prop_kg + cp_mass_kg
    dv_ks = stage_delta_v(ks_isp_s, ks_struct_kg, ks_prop_kg, payload_for_ks)

    # DSPS carries only CP
    payload_for_dsps = cp_mass_kg
    dv_dsps = stage_delta_v(dsps_isp_s, dsps_struct_kg, dsps_prop_kg, payload_for_dsps)

    return dv_ks + dv_dsps


# =========================================
# 5. CP_max 계산 함수 (Isp, prop mass를 바꿔가며)
# =========================================

def cp_max_for_config(
    mode: str,
    margin_m_per_s: float,
    ks_isp_s: float,
    dsps_isp_s: float,
    ks_struct_kg: float,
    dsps_struct_kg: float,
    ks_prop_kg: float,
    dsps_prop_kg: float,
    nuri_leo_limit_kg: float = 3200.0,
    cp_search_max_kg: int = 800,
    cp_step_kg: int = 1,
    leo_alt_km: float = 300.0,
    venus_orbit_alt_km: float = 300.0,
):
    """
    주어진 Isp/prop 구조에서 CP_max를 찾는다.

    - direct_entry  : TVI + margin만 요구
    - orbit_insertion: TVI + VOI + margin까지 요구

    반환:
        (cp_max_kg, dv_req_m_s, note_str)
        cp_max_kg == -1 이면, CP=0kg도 ΔV 요구를 만족 못한다는 뜻.
    """
    dv_info = dv_required(
        mode=mode,
        leo_alt_km=leo_alt_km,
        venus_orbit_alt_km=venus_orbit_alt_km,
        margin_m_per_s=margin_m_per_s,
    )
    dv_req = dv_info["dv_total"]

    # 먼저 CP=0 kg일 때도 안 되는지 확인
    m_leo_0 = ks_struct_kg + ks_prop_kg + dsps_struct_kg + dsps_prop_kg + 0.0
    dv_0 = total_dv_for_cp(
        0.0,
        ks_isp_s,
        ks_struct_kg,
        ks_prop_kg,
        dsps_isp_s,
        dsps_struct_kg,
        dsps_prop_kg,
    )

    if dv_0 < dv_req:
        # ΔV가 애초에 부족
        return -1.0, dv_req, "ΔV 부족: CP=0 kg도 요구 ΔV를 만족 못함."

    # CP를 늘려가며 DV + LEO mass 두 조건을 동시에 만족하는 최대값 탐색
    cp_max = 0.0
    for cp in range(0, cp_search_max_kg + 1, cp_step_kg):
        # LEO mass check
        m_leo = ks_struct_kg + ks_prop_kg + dsps_struct_kg + dsps_prop_kg + cp
        if m_leo > nuri_leo_limit_kg:
            # 누리 LEO 한계 초과 → 더 큰 CP는 볼 필요 없음 (mass 제한)
            break

        dv_total = total_dv_for_cp(
            cp,
            ks_isp_s,
            ks_struct_kg,
            ks_prop_kg,
            dsps_isp_s,
            dsps_struct_kg,
            dsps_prop_kg,
        )
        if dv_total >= dv_req:
            cp_max = cp
        else:
            # CP 증가하면 ΔV 더 줄어드니까, 여기서 멈춰도 됨
            break

    note = "성공: 이 CP까지는 ΔV+LEO 한계 모두 만족."
    return cp_max, dv_req, note


# =========================================
# 6. Isp & prop mass trade-off 스터디
# =========================================

def trade_isp_and_prop_for_cp_max():
    # ---- 기본 KARI 구조질량 ----
    KS_STRUCT_BASE = 450.0   # kg
    DSPS_STRUCT_BASE = 180.0 # kg

    # ---- 기본 KARI 추진제 질량 ----
    KS_PROP_BASE = 2050.0    # kg
    DSPS_PROP_BASE = 120.0   # kg

    # ---- 누리호 LEO 한계 ----
    NURI_LEO_LIMIT = 3200.0  # kg

    # ---- 설정 ----
    mode = "orbit_insertion"     # 또는 "orbit_insertion"
    margin = 300.0            # [m/s]
    leo_alt_km = 300.0
    venus_alt_km = 300.0

    # 스윕할 Isp 및 prop scale 리스트
    ks_isp_list = [340.0, 360.0, 380.0]     # KS Isp 후보
    dsps_isp_list = [280.0, 290.0, 300.0]   # DSPS Isp 후보
    prop_scale_list = [0.8, 1.0, 1.2]       # 기본 추진제 질량 대비 배율

    print("===== Isp & prop mass trade study for CP_max =====")
    print(f"mode = {mode}, margin = {margin} m/s\n")
    header = "KS_Isp | DSPS_Isp | prop_scale | CP_max [kg] | note"
    print(header)

    for ks_isp in ks_isp_list:
        for dsps_isp in dsps_isp_list:
            for prop_scale in prop_scale_list:
                ks_prop = KS_PROP_BASE * prop_scale
                dsps_prop = DSPS_PROP_BASE * prop_scale

                cp_max, dv_req, note = cp_max_for_config(
                    mode=mode,
                    margin_m_per_s=margin,
                    ks_isp_s=ks_isp,
                    dsps_isp_s=dsps_isp,
                    ks_struct_kg=KS_STRUCT_BASE,
                    dsps_struct_kg=DSPS_STRUCT_BASE,
                    ks_prop_kg=ks_prop,
                    dsps_prop_kg=dsps_prop,
                    nuri_leo_limit_kg=NURI_LEO_LIMIT,
                    cp_search_max_kg=800,
                    cp_step_kg=1,
                    leo_alt_km=leo_alt_km,
                    venus_orbit_alt_km=venus_alt_km,
                )

                print(f"{ks_isp:6.1f} | {dsps_isp:9.1f} |"
                      f" {prop_scale:9.2f} | {cp_max:11.1f} | {note}")

    print("\n(참고) dv_req (mode 기준)은 위에서 cp_max_for_config 안에서 계산됨.")


if __name__ == "__main__":
    trade_isp_and_prop_for_cp_max()
