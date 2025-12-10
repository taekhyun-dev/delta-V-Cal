import numpy as np

# 1. Mission Parameters
mass_arrival = 700      # kg (Delivered by KS)
isp_dsps = 290          # sec (Propulsion efficiency)
g0 = 9.81
mass_dsps_dry = 180     # kg (Bus/Structure mass included in Main S/C)

# 2. Orbit Insertion Delta V (Elliptical Equatorial)
# Conservative estimate for VOI + Phasing into formation
delta_v_target = 1400   # m/s

# 3. Calculate Total Mass in Orbit (M_final)
# Rocket Equation: M_final = M_initial / exp(dV / (Isp*g0))
mass_ratio = np.exp(delta_v_target / (isp_dsps * g0))
mass_orbit_total = mass_arrival / mass_ratio
mass_fuel_burned = mass_arrival - mass_orbit_total

# 4. Scenario Analysis: Main S/C Mass Calculation
# Scenarios for Sub-S/C mass (each)
sub_sc_masses = [10, 20, 40] # kg (6U, 12U, Micro)

print(f"--- Venus Formation Flying Mission Design ---")
print(f"Target Orbit: Elliptical Equatorial (Delta V ~ {delta_v_target} m/s)")
print(f"Arrival Mass: {mass_arrival} kg")
print(f"Fuel Consumed for VOI: {mass_fuel_burned:.2f} kg")
print(f"Total Mass Available in Orbit: {mass_orbit_total:.2f} kg")
print("-" * 50)
print(f"{'Sub-S/C (each)':<15} | {'Sub Total (x2)':<15} | {'Main S/C Total':<15} | {'Main S/C Payload*'}")
print("-" * 50)

for sub_mass in sub_sc_masses:
    sub_total = sub_mass * 2
    main_sc_total = mass_orbit_total - sub_total
    
    # Main S/C Payload = Main Total - DSPS Structure (180kg)
    main_payload = main_sc_total - mass_dsps_dry
    
    print(f"{sub_mass} kg {'':<11} | {sub_total} kg {'':<11} | {main_sc_total:.2f} kg {'':<8} | {main_payload:.2f} kg")

print("-" * 50)
print("* 'Main S/C Total' includes the DSPS bus structure (180kg).")
print("* 'Main S/C Payload' is pure capacity for Science Instruments + Avionics.")