import math

# ============================
# 1. 상수 정의
# ============================

G0 = 9.80665  # m/s^2

MU_SUN = 1.32712440018e20   # m^3/s^2
MU_EARTH = 3.986004418e14   # m^3/s^2
R_EARTH = 6371e3            # m
AU = 1.495978707e11         # m

R_VENUS = 6051.8e3          # m
MU_VENUS = 3.24859e14       # m^3/s^2

R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# ============================
# 2. 기본 궤도역학 함수
# ============================

def earth_venus_hohmann():
    """
    태양 중심 Hohmann 전이 (지구 1 AU -> 금성 0.723 AU) 파라미터 계산.
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
    LEO(원궤도) -> 주어진 v_inf를 가진 탈출 궤도까지 ΔV 계산.
    """
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv, v_circ, v_esc, v_peri_hyp


def dv_venus_capture(v_inf_venus, venus_orbit_alt_km=300.0):
    """
    금성 도착 시 v_inf_venus를 가진 하이퍼볼릭 궤도에서
    금성 고도 venus_orbit_alt_km 원궤도로 캡처하는 ΔV 계산.
    """
    r_p = R_VENUS + venus_orbit_alt_km * 1e3

    v_circ = math.sqrt(MU_VENUS / r_p)
    v_esc = math.sqrt(2.0 * MU_VENUS / r_p)
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)

    dv_capture = v_peri_hyp - v_circ
    return dv_capture, v_circ, v_esc, v_peri_hyp


# ============================
# 3. ΔV 요구량 계산
# ============================

def venus_dv_required(
    mode="direct_entry",        # "direct_entry" 또는 "orbit_insertion"
    leo_alt_km=300.0,
    venus_orbit_alt_km=300.0,
    margin_m_per_s=300.0
):
    """
    mode에 따라 요구 ΔV 계산:
    - direct_entry: TVI + margin (금성 대기 진입, VOI 없음)
    - orbit_insertion: TVI + VOI + margin (금성 궤도 안착)
    """
    hv = earth_venus_hohmann()

    # TVI (LEO -> 전이궤도)
    dv_tvi, v_circ_leo, v_esc_leo, v_peri_hyp_leo = dv_leo_to_vinf(
        hv["v_inf_earth"], leo_alt_km
    )

    dv_voi = 0.0
    if mode == "orbit_insertion":
        dv_voi, v_circ_v, v_esc_v, v_peri_hyp_v = dv_venus_capture(
            hv["v_inf_venus"], venus_orbit_alt_km
        )

    dv_total = dv_tvi + dv_voi + margin_m_per_s

    return {
        "mode": mode,
        "dv_tvi": dv_tvi,
        "dv_voi": dv_voi,
        "dv_margin": margin_m_per_s,
        "dv_total": dv_total,
    }


# ============================
# 4. KS / DSPS ΔV 계산
# ============================

def stage_delta_v(isp_s, m_struct, m_prop, m_payload):
    """
    단일 스테이지 ΔV:
    isp_s    : Isp [s]
    m_struct : 구조체 질량 [kg]
    m_prop   : 추진제 질량 [kg]
    m_payload: 이 스테이지가 들고 있는 상단+캡슐 질량 [kg]
    """
    if m_prop <= 0.0:
        return 0.0

    m0 = m_struct + m_prop + m_payload
    mf = m_struct + m_payload

    if m0 <= mf:
        return 0.0

    vex = isp_s * G0
    return vex * math.log(m0 / mf)


def total_dv_for_cp(cp_mass_kg, ks_cfg, dsps_cfg):
    """
    주어진 KS/DSPS 설계에서 CP 질량 cp_mass_kg일 때
    KS + DSPS가 제공하는 총 ΔV 계산.
    ks_cfg, dsps_cfg는 dict로:
        {"isp": .., "m_struct": .., "m_prop": ..}
    """
    # KS가 들고 가는 payload = DSPS 전체 + CP
    payload_for_ks = (
        dsps_cfg["m_struct"] +
        dsps_cfg["m_prop"] +
        cp_mass_kg
    )
    dv_ks = stage_delta_v(
        ks_cfg["isp"],
        ks_cfg["m_struct"],
        ks_cfg["m_prop"],
        payload_for_ks
    )

    # DSPS는 CP만 들고 burn
    dv_dsps = stage_delta_v(
        dsps_cfg["isp"],
        dsps_cfg["m_struct"],
        dsps_cfg["m_prop"],
        cp_mass_kg
    )

    return dv_ks + dv_dsps


# ============================
# 5. 특정 설계에서 CP 최대값 찾기
# ============================

def find_cp_max_for_config(
    dv_required_m_per_s,
    ks_cfg,
    dsps_cfg,
    nuri_leo_limit_kg=3200.0,
    cp_search_max_kg=800,
    cp_step_kg=1
):
    """
    주어진 KS/DSPS 설계 + ΔV 요구량에서,
    - 누리 LEO 질량 한계
    - ΔV 요구량 충족
    를 동시에 만족하는 CP 최대 질량을 찾는다.

    반환: (cp_max, reason)
      - cp_max < 0 이면 완전 불가능(0 kg도 ΔV 미달)
    """
    cp_max = -1.0

    # KS, DSPS total mass (CP 제외)
    m_ks_total = ks_cfg["m_struct"] + ks_cfg["m_prop"]
    m_dsps_total = dsps_cfg["m_struct"] + dsps_cfg["m_prop"]

    for cp in range(0, int(cp_search_max_kg) + 1, cp_step_kg):
        m_leo = m_ks_total + m_dsps_total + cp
        if m_leo > nuri_leo_limit_kg:
            continue  # 누리 한계 초과

        dv_total = total_dv_for_cp(cp, ks_cfg, dsps_cfg)
        if dv_total >= dv_required_m_per_s:
            cp_max = cp  # 조건 만족, 더 큰 cp를 계속 탐색

    if cp_max < 0:
        reason = "ΔV 부족: CP=0 kg도 요구 ΔV를 만족 못함."
    else:
        reason = "성공: 이 CP까지는 ΔV+LEO 한계 모두 만족."

    return cp_max, reason


# ============================
# 6. 트레이드 스터디
# ============================

def trade_study():
    # --- KARI baseline KS/DSPS ---
    ks_baseline = {
        "isp": 360.0,
        "m_struct": 450.0,
        "m_prop": 2050.0,
    }
    dsps_baseline = {
        "isp": 290.0,
        "m_struct": 180.0,
        "m_prop": 120.0,
    }

    nuri_leo_limit = 3200.0

    print("===== [A] margin 변화에 따른 CP_max (direct_entry / orbit_insertion) =====")
    for mode in ["direct_entry", "orbit_insertion"]:
        print(f"\n--- mode = {mode} ---")
        print("margin [m/s] | dv_req [km/s] | CP_max [kg] | note")
        for margin in [0.0, 100.0, 200.0, 300.0, 400.0]:
            dv_req = venus_dv_required(
                mode=mode,
                leo_alt_km=300.0,
                venus_orbit_alt_km=300.0,
                margin_m_per_s=margin
            )
            cp_max, reason = find_cp_max_for_config(
                dv_required_m_per_s=dv_req["dv_total"],
                ks_cfg=ks_baseline,
                dsps_cfg=dsps_baseline,
                nuri_leo_limit_kg=nuri_leo_limit,
                cp_search_max_kg=800,
                cp_step_kg=1
            )
            print(f"{margin:10.1f} | {dv_req['dv_total']/1000:11.3f} | {cp_max:10.1f} | {reason}")

    print("\n\n===== [B] KS/DSPS 구조 변경에 따른 CP_max (direct_entry, margin=300 m/s) =====")
    mode = "direct_entry"
    dv_req = venus_dv_required(
        mode=mode,
        leo_alt_km=300.0,
        venus_orbit_alt_km=300.0,
        margin_m_per_s=300.0
    )
    dv_req_val = dv_req["dv_total"]

    # 예시: KS Isp, DSPS Isp, 추진제 스케일을 바꿔보는 간단한 grid
    ks_isp_list = [340.0, 360.0, 380.0]     # KS Isp 변화
    dsps_isp_list = [280.0, 290.0, 300.0]   # DSPS Isp 변화
    prop_scale_list = [0.8, 1.0, 1.2]       # 추진제 질량 스케일

    print("KS_Isp | DSPS_Isp | prop_scale | CP_max [kg] | note")
    for ks_isp in ks_isp_list:
        for dsps_isp in dsps_isp_list:
            for scale in prop_scale_list:
                ks_cfg = {
                    "isp": ks_isp,
                    "m_struct": ks_baseline["m_struct"],
                    "m_prop": ks_baseline["m_prop"] * scale,
                }
                dsps_cfg = {
                    "isp": dsps_isp,
                    "m_struct": dsps_baseline["m_struct"],
                    "m_prop": dsps_baseline["m_prop"] * scale,
                }

                cp_max, reason = find_cp_max_for_config(
                    dv_required_m_per_s=dv_req_val,
                    ks_cfg=ks_cfg,
                    dsps_cfg=dsps_cfg,
                    nuri_leo_limit_kg=nuri_leo_limit,
                    cp_search_max_kg=800,
                    cp_step_kg=1
                )

                print(f"{ks_isp:6.1f} | {dsps_isp:9.1f} | {scale:10.2f} | {cp_max:10.1f} | {reason}")


if __name__ == "__main__":
    trade_study()
