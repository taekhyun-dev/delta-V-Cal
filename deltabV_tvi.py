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
R_VENUS_ORBIT = 0.723 * AU  # 단순 원궤도 가정


# ============================
# 2. Hohmann Earth→Venus transfer
# ============================

def earth_venus_hohmann():
    """
    태양 중심 Hohmann 전이 (지구 1 AU -> 금성 0.723 AU)
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


def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    """
    LEO(원궤도)에서 주어진 v_inf를 갖는 탈출 궤도로 가는 ΔV
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv, v_circ, v_esc, v_peri_hyp


def venus_direct_entry_dv_required(leo_alt_km=300.0, margin_m_per_s=300.0):
    """
    금성 '대기 진입'만 고려한 ΔV 요구:
    - LEO -> 금성 전이(TVI)
    - + margin (mid-course, targeting 등)
    - 금성 궤도 포획(VOI) 없음, 지구 귀환 없음
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
    단일 스테이지 로켓 방정식 ΔV
    isp_s    : Isp [s]
    m_struct : 구조 질량 [kg]
    m_prop   : 추진제 질량 [kg]
    m_payload: 이 스테이지가 들고 있는 페이로드 질량 [kg]
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
    주어진 CP 질량(cp_mass_kg)에 대해
    Kick Stage + DSPS가 제공할 수 있는 총 ΔV [m/s]
    (KARI Table 2 기반 스펙 사용)
    """
    # ---- Kick Stage (KS) ----
    ISP_KS = 360.0
    M_KS_STRUCT = 450.0
    M_KS_PROP = 2050.0

    # ---- Deep Space Propulsion System (DSPS) ----
    ISP_DSPS = 290.0
    M_DSPS_STRUCT = 180.0
    M_DSPS_PROP = 120.0

    # KS: DSPS(구조+추진제) + CP를 페이로드로 들고 TVI 수행
    payload_for_ks = M_DSPS_STRUCT + M_DSPS_PROP + cp_mass_kg
    dv_ks = stage_delta_v(ISP_KS, M_KS_STRUCT, M_KS_PROP, payload_for_ks)

    # DSPS: CP만 들고 mid-course/entry targeting ΔV 수행
    payload_for_dsps = cp_mass_kg
    dv_dsps = stage_delta_v(ISP_DSPS, M_DSPS_STRUCT, M_DSPS_PROP, payload_for_dsps)

    return dv_ks + dv_dsps


# ============================
# 4. Nuri + KS + DSPS + Capsule (direct entry case)
# ============================

def nuri_venus_direct_entry_scenario():
    # ---- 환경 / 미션 파라미터 ----
    LEO_ALT_KM = 300.0
    DV_MARGIN = 300.0  # m/s

    # 누리호 LEO 성능 (KARI 시나리오)
    M_NURI_LEO_LIMIT = 3200.0  # kg

    # KS / DSPS 질량 (KARI Table 2)
    M_KS_TOTAL = 2500.0  # struct+prop
    M_DSPS_TOTAL = 300.0

    # 1) ΔV 요구 (TVI + margin)
    dv_req = venus_direct_entry_dv_required(
        leo_alt_km=LEO_ALT_KM,
        margin_m_per_s=DV_MARGIN
    )
    dv_tvi = dv_req["dv_tvi"]
    dv_total_req = dv_req["dv_total"]

    print("===== Venus DIRECT ENTRY ΔV budget =====")
    print(f"LEO altitude          : {LEO_ALT_KM:.1f} km")
    print(f"ΔV (TVI, Earth dep.)  : {dv_tvi/1000:.3f} km/s")
    print(f"ΔV margin             : {DV_MARGIN/1000:.3f} km/s")
    print(f"TOTAL ΔV required     : {dv_total_req/1000:.3f} km/s\n")

    # 2) KARI baseline CP = 400 kg일 때 ΔV 확인
    baseline_cp = 400.0
    dv_baseline = total_dv(baseline_cp)
    print("===== Baseline (KARI-like) performance =====")
    print(f"Baseline CP mass      : {baseline_cp:.1f} kg")
    print(f"Total ΔV (KS+DSPS)    : {dv_baseline/1000:.3f} km/s\n")

    # 3) CP 최대 질량 스캔
    cp_max = 0.0
    for cp in range(0, 801):  # 0~800 kg 범위에서 1 kg 단위로 탐색
        m_leo = M_KS_TOTAL + M_DSPS_TOTAL + cp
        if m_leo > M_NURI_LEO_LIMIT:
            continue  # 누리 LEO 한계 초과 → 불가능

        dv_total = total_dv(cp)
        if dv_total >= dv_total_req:  # 요구 ΔV 만족?
            cp_max = cp

    print("===== CP mass search (direct entry only) =====")
    if cp_max <= 0.0:
        print("-> Even with 0 kg capsule, required ΔV cannot be met.")
    else:
        print(f"-> Max CP mass satisfying ΔV + LEO constraints: {cp_max:.1f} kg")
        dv_at_cp_max = total_dv(cp_max)
        print(f"   ΔV at that CP mass                        : {dv_at_cp_max/1000:.3f} km/s")
        print(f"   ΔV margin                                 : "
              f"{(dv_at_cp_max - dv_total_req):.1f} m/s")


if __name__ == '__main__':
    nuri_venus_direct_entry_scenario()
