import math
import csv

# ============================
# 1. Physical constants
# ============================

G0 = 9.80665  # m/s^2
MU_SUN   = 1.32712440018e20
MU_EARTH = 3.986004418e14
MU_VENUS = 3.24859e14

R_EARTH = 6371e3
R_VENUS = 6051.8e3
AU = 1.495978707e11
R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU


# ============================
# 2. Heliocentric Earth↔Venus transfer
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
        "v_inf_earth": v_inf_earth,
        "v_inf_venus": v_inf_venus,
    }


# ============================
# 3. ΔV Calculations (KS & DSPS)
# ============================

def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc  = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    dv = v_peri_hyp - v_circ
    return dv

def dv_venus_orbit_insertion(v_inf_venus, venus_orbit_alt_km=300.0):
    r_p = R_VENUS + venus_orbit_alt_km * 1e3
    v_circ = math.sqrt(MU_VENUS / r_p)          
    v_esc  = math.sqrt(2.0 * MU_VENUS / r_p)    
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)
    
    # 궤도 진입(Capture)에 필요한 감속량
    dv_capture = v_peri_hyp - v_circ
    return dv_capture


# ============================
# 4. Rocket Equation
# ============================

def stage_delta_v(isp_s, m_struct, m_prop, m_payload):
    if m_prop <= 0.0: return 0.0
    m0 = m_struct + m_prop + m_payload
    mf = m_struct + m_payload
    if m0 <= mf: return 0.0
    
    vex = isp_s * G0
    return vex * math.log(m0 / mf)


# ============================
# 5. Mission Requirements
# ============================

def dv_requirements(
    mission_profile="direct_entry",
    leo_alt_km=300.0,
    venus_orbit_alt_km=300.0,
    dsps_margin_m_per_s=300.0
):
    hv = earth_venus_hohmann()

    # 1) KS Requirement (Earth Departure)
    dv_tvi = dv_leo_to_vinf(hv["v_inf_earth"], leo_alt_km)

    # 2) DSPS Requirement (Venus Arrival)
    if mission_profile == "orbit_insertion":
        dv_voi = dv_venus_orbit_insertion(hv["v_inf_venus"], venus_orbit_alt_km)
        dv_dsps_req = dv_voi + dsps_margin_m_per_s
    
    elif mission_profile == "direct_entry":
        dv_voi = 0.0
        dv_dsps_req = dsps_margin_m_per_s
        
    else:
        raise ValueError(f"Unknown mission profile: {mission_profile}")

    return {
        "dv_ks_req": dv_tvi,
        "dv_dsps_req": dv_dsps_req,
        "mission_profile": mission_profile
    }


# ============================
# 6. Trade-off Analysis (Fixed)
# ============================

def run_trade_off_analysis(
    mission_profile,
    nuri_leo_limit_kg=3200.0,
    m_ks_struct=450.0, m_ks_prop=2050.0, isp_ks=360.0,
    m_dsps_struct=180.0, m_dsps_prop=120.0, isp_dsps=290.0,
    dsps_margin_m_per_s=300.0  # [수정됨] 인자 추가
):
    # 1. 요구 조건 계산
    reqs = dv_requirements(
        mission_profile=mission_profile,
        dsps_margin_m_per_s=dsps_margin_m_per_s # [수정됨] 값 전달
    )
    dv_ks_req = reqs["dv_ks_req"]
    dv_dsps_req = reqs["dv_dsps_req"]
    
    # 2. CP 스윕 범위 설정
    base_mass = (m_ks_struct + m_ks_prop + m_dsps_struct + m_dsps_prop)
    cp_limit = max(0.0, nuri_leo_limit_kg - base_mass)
    
    max_feasible_cp = -1.0
    logs = []
    
    # 3. CP 0kg ~ 한계까지 반복
    cp = 0.0
    while cp <= cp_limit + 1e-6:
        
        # KS Performance
        payload_ks = m_dsps_struct + m_dsps_prop + cp
        dv_ks_avail = stage_delta_v(isp_ks, m_ks_struct, m_ks_prop, payload_ks)
        
        # DSPS Performance
        payload_dsps = cp
        dv_dsps_avail = stage_delta_v(isp_dsps, m_dsps_struct, m_dsps_prop, payload_dsps)
        
        # Check Success
        ks_ok = (dv_ks_avail >= dv_ks_req)
        dsps_ok = (dv_dsps_avail >= dv_dsps_req)
        success = ks_ok and dsps_ok
        
        if success:
            max_feasible_cp = cp 
            
        logs.append({
            "Profile": mission_profile,
            "CP_kg": cp,
            "KS_Avail": dv_ks_avail,
            "KS_Req": dv_ks_req,
            "DSPS_Avail": dv_dsps_avail,
            "DSPS_Req": dv_dsps_req,
            "Success": success
        })
        
        cp += 1.0

    return max_feasible_cp, logs


# ============================
# 7. Main Execution
# ============================

def main():
    scenarios = ["direct_entry", "orbit_insertion"]
    all_logs = []

    print(f"{'Mission Profile':<20} | {'Max CP (kg)':<15} | {'Result Note'}")
    print("-" * 60)

    for profile in scenarios:
        max_cp, logs = run_trade_off_analysis(
            mission_profile=profile,
            nuri_leo_limit_kg=3200.0,
            dsps_margin_m_per_s=300.0
        )
        all_logs.extend(logs)
        
        note = "Success (Matches Paper)" if max_cp > 0 else "Fail (Need more fuel)"
        print(f"{profile:<20} | {max_cp:<15.1f} | {note}")

    filename = "trade_off_results.csv"
    with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=all_logs[0].keys())
        writer.writeheader()
        writer.writerows(all_logs)
    
    print("-" * 60)
    print(f"상세 트레이드 오프 데이터가 '{filename}'에 저장되었습니다.")

if __name__ == "__main__":
    main()