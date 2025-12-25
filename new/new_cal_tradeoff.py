import math
import csv

# ============================
# 1. Physical constants
# ============================

G0 = 9.80665
MU_SUN   = 1.32712440018e20
MU_EARTH = 3.986004418e14
MU_VENUS = 3.24859e14

R_EARTH = 6371e3
R_VENUS = 6051.8e3
AU = 1.495978707e11
R_EARTH_ORBIT = AU
R_VENUS_ORBIT = 0.723 * AU

# ============================
# 2. Heliocentric Calculation
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

    return {"v_inf_earth": v_inf_earth, "v_inf_venus": v_inf_venus}

# ============================
# 3. Delta-V Logic
# ============================

def dv_leo_to_vinf(v_inf, leo_alt_km=300.0):
    r_leo = R_EARTH + leo_alt_km * 1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc  = math.sqrt(2.0 * MU_EARTH / r_leo)
    v_peri_hyp = math.sqrt(v_inf**2 + v_esc**2)
    return v_peri_hyp - v_circ

def dv_venus_orbit_insertion(v_inf_venus, venus_orbit_alt_km=300.0):
    r_p = R_VENUS + venus_orbit_alt_km * 1e3
    v_circ = math.sqrt(MU_VENUS / r_p)
    v_esc  = math.sqrt(2.0 * MU_VENUS / r_p)
    v_peri_hyp = math.sqrt(v_inf_venus**2 + v_esc**2)
    return v_peri_hyp - v_circ

def stage_delta_v(isp_s, m_struct, m_prop, m_payload):
    if m_prop <= 0.0: return 0.0
    m0 = m_struct + m_prop + m_payload
    mf = m_struct + m_payload
    if m0 <= mf: return 0.0
    vex = isp_s * G0
    return vex * math.log(m0 / mf)

# ============================
# 4. Mission Requirements (Updated)
# ============================

def get_dv_requirements(mission_profile, leo_alt_km, venus_orbit_alt_km, dsps_margin_m_per_s):
    hv = earth_venus_hohmann()
    dv_tvi = dv_leo_to_vinf(hv["v_inf_earth"], leo_alt_km)

    # --- Profile별 요구조건 분기 ---
    if mission_profile == "orbit_insertion":
        # 1. 순수 추진 궤도 진입: (감속 dV) + 마진
        dv_voi = dv_venus_orbit_insertion(hv["v_inf_venus"], venus_orbit_alt_km)
        dv_dsps_req = dv_voi + dsps_margin_m_per_s
        
    elif mission_profile == "aerocapture":
        # 2. 공력 포획: (근지점 상승 기동 PRM) + 마진
        # PRM은 통상 100~200 m/s 소요. 보수적으로 200 m/s 적용
        dv_prm = 200.0 
        dv_dsps_req = dv_prm + dsps_margin_m_per_s
        
    elif mission_profile == "direct_entry":
        # 3. 대기 직접 진입: 0 + 마진 (순수 TCM)
        dv_dsps_req = dsps_margin_m_per_s
        
    else:
        raise ValueError("Unknown profile")

    return dv_tvi, dv_dsps_req

# ============================
# 5. Trade-off Optimization
# ============================

def optimize_mass_distribution(
    mission_profile,
    nuri_leo_limit_kg=3200.0,
    m_ks_struct=450.0, m_ks_prop=2050.0, isp_ks=360.0,
    m_dsps_struct=180.0, isp_dsps=290.0,
    dsps_margin_m_per_s=300.0
):
    dv_ks_req, dv_dsps_req = get_dv_requirements(
        mission_profile, 300.0, 300.0, dsps_margin_m_per_s
    )

    m_ks_wet = m_ks_struct + m_ks_prop
    available_mass_for_upper = nuri_leo_limit_kg - m_ks_wet
    max_possible_prop = available_mass_for_upper - m_dsps_struct

    best_result = {
        "cp_kg": -1.0, "dsps_prop_kg": -1.0, 
        "note": "Fail", "dv_dsps_cap": 0.0
    }
    
    # 추진제 양을 0부터 늘려가며 탐색 (Min Propellant -> Max CP)
    for m_prop in range(0, int(max_possible_prop) + 1):
        cp_current = available_mass_for_upper - m_dsps_struct - m_prop
        if cp_current < 0: break

        # KS Check
        payload_for_ks = m_dsps_struct + m_prop + cp_current
        dv_ks_avail = stage_delta_v(isp_ks, m_ks_struct, m_ks_prop, payload_for_ks)
        
        # DSPS Check
        payload_for_dsps = cp_current
        dv_dsps_avail = stage_delta_v(isp_dsps, m_dsps_struct, m_prop, payload_for_dsps)
        
        if dv_ks_avail >= dv_ks_req and dv_dsps_avail >= dv_dsps_req:
            best_result = {
                "cp_kg": cp_current,
                "dsps_prop_kg": m_prop,
                "note": "Success",
                "dv_dsps_cap": dv_dsps_avail
            }
            break # 찾았으면 루프 종료 (최대 CP)

    return best_result, dv_dsps_req

# ============================
# 6. Main Execution
# ============================

def main():
    margins = [0, 50, 100, 150, 200, 250, 300]
    profiles = ["direct_entry", "aerocapture", "orbit_insertion"] # 3개 프로파일
    
    all_results = []

    print(f"{'Margin':<6} | {'Profile':<18} | {'Max CP (kg)':<12} | {'Req Prop':<10} | {'Req dV'}")
    print("-" * 75)

    for margin in margins:
        for profile in profiles:
            res, req_dv = optimize_mass_distribution(
                mission_profile=profile,
                nuri_leo_limit_kg=3200.0,
                dsps_margin_m_per_s=float(margin)
            )
            
            all_results.append({
                "Margin": margin,
                "Profile": profile,
                "Max_CP_kg": res["cp_kg"],
                "Req_Prop_kg": res["dsps_prop_kg"],
                "Req_dV": req_dv,
                "Note": res["note"]
            })
            
            print(f"{margin:<6} | {profile:<18} | {res['cp_kg']:<12.1f} | {res['dsps_prop_kg']:<10.1f} | {req_dv:.1f}")

    # CSV 저장
    filename = "three_profile_tradeoff.csv"
    with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
        fieldnames = ["Margin", "Profile", "Max_CP_kg", "Req_Prop_kg", "Req_dV", "Note"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print("-" * 75)
    print(f"3가지 프로파일 분석 결과가 '{filename}'에 저장되었습니다.")

if __name__ == "__main__":
    main()