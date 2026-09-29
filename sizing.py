from math import pi, sqrt

MOMENT_HIGH_SPEED = 9_000  # Nm, hinge moment on one rudder surface
RATE_HIGH_SPEED = 0.30  # rad/s
DEFLECTION_HIGH_SPEED = 10  # deg

MOMENT_LOW_SPEED = 1_000  # Nm
RATE_LOW_SPEED = 0.20  # rad/s
DEFLECTION_LOW_SPEED = 60  # deg

MECHANISMS = 2  # actuators per surface, sharing the moment
SURFACES = 4  # upper + lower, left + right

# lever arm is a jacobian...
LEVER_ARM_HIGH_SPEED = 78.8e-3  # m/rad, at the high-speed deflection
LEVER_ARM_LOW_SPEED = 54.7e-3  # m/rad, at the low-speed deflection
LEVER_ARM_PEAK = 78.8e-3  # m/rad, the largest anywhere in the travel
STROKE = 90.5e-3  # m, rod travel over the full range
HORN_LENGTH = 50e-3  # m, hinge to the pushrod pin
LOAD_PATHS = 1  # ball links sharing the load at the horn

SUPPLY_PRESSURE = 28e6  # Pa
ROD_RATIO = 0.20  # rod diameter / bore diameter
EFFICIENCY = 0.85  # seal friction and joint losses
BULK_MODULUS = 1.5e9  # Pa


def rod_force(moment, lever_arm):
    """Force one actuator must produce to hold this moment.

    The moment is shared between the mechanisms, and the lever arm converts
    the remaining moment into a force.
    """
    moment_per_mechanism = moment / MECHANISMS
    return moment_per_mechanism / lever_arm


def piston_areas(force):
    """(cap side, rod side) piston areas needed to make this force.

    Retracting is the weak stroke, because the rod eats into that face, so
    the cylinder is sized on the annulus.
    """
    annulus_area = force / (EFFICIENCY * SUPPLY_PRESSURE)
    bore_area = annulus_area / (1 - ROD_RATIO**2)
    return bore_area, annulus_area


def bore_diameter(bore_area):
    """Cylinder bore diameter for this piston area."""
    return 2 * sqrt(bore_area / pi)


def rod_speed(lever_arm, rate):
    """How fast the rod moves when the rudder turns at this rate."""
    return lever_arm * rate


def flow_litres_per_minute(area, speed):
    """Oil displaced by a piston of this area moving at this speed."""
    return area * speed * 60_000


def pushrod_load(moment):
    """TOTAL force in the link between the horn and the bellcrank."""
    moment_per_mechanism = moment / MECHANISMS
    return moment_per_mechanism / HORN_LENGTH


def load_per_joint(total_load):
    """What one ball link carries, if the load is split between several."""
    return total_load / LOAD_PATHS


def oil_stiffness_at_rod(bore_area, annulus_area):
    """Spring rate of the trapped oil, felt at the rod."""
    column_length = STROKE / 2
    return BULK_MODULUS * (bore_area + annulus_area) / column_length


def stiffness_at_hinge(rod_stiffness, lever_arm):
    """The same spring, felt at the rudder hinge instead of at the rod."""
    return MECHANISMS * lever_arm**2 * rod_stiffness


def main():
    # what each flight case demands
    force_high = rod_force(MOMENT_HIGH_SPEED, LEVER_ARM_HIGH_SPEED)
    force_low = rod_force(MOMENT_LOW_SPEED, LEVER_ARM_LOW_SPEED)

    speed_high = rod_speed(LEVER_ARM_PEAK, RATE_HIGH_SPEED)
    speed_low = rod_speed(LEVER_ARM_PEAK, RATE_LOW_SPEED)

    # size the cylinder on the worst of them
    design_force = max(force_high, force_low)
    design_speed = max(speed_high, speed_low)

    bore_area, annulus_area = piston_areas(design_force)
    bore = bore_diameter(bore_area)

    print("FLOW AND POWER")

    flow_cap = flow_litres_per_minute(bore_area, design_speed)
    flow_rod = flow_litres_per_minute(annulus_area, design_speed)
    flow_surface = flow_cap * MECHANISMS
    flow_system = flow_surface * SURFACES
    power_system = SUPPLY_PRESSURE * flow_system / 60_000

    print(f"  per actuator, cap side   {flow_cap:8.2f} L/min")
    print(f"  per actuator, rod side   {flow_rod:8.2f} L/min")
    print(f"  per surface              {flow_surface:8.2f} L/min")
    print(f"  whole system             {flow_system:8.2f} L/min")
    print(f"  hydraulic power          {power_system / 1000:8.2f} kW")
    print()

    print()
    print("PUSHROD LOAD")

    load = pushrod_load(MOMENT_HIGH_SPEED)
    per_joint = load_per_joint(load)

    print(f"  horn length              {HORN_LENGTH * 1000:8.1f} mm")
    print(f"  pushrod load, total      {load / 1000:8.1f} kN")
    print(f"  load paths               {LOAD_PATHS:8d}")
    print(f"  per ball link            {per_joint / 1000:8.1f} kN")
    print(f"  actuator rod load        {design_force / 1000:8.1f} kN")
    print(f"  ratio, pushrod / rod     {load / design_force:8.1f} x")
    print()

    print()
    print("STIFFNESS")

    rod_stiffness = oil_stiffness_at_rod(bore_area, annulus_area)
    hinge_stiffness = stiffness_at_hinge(rod_stiffness, LEVER_ARM_PEAK)

    print(f"  oil column at the rod    {rod_stiffness / 1e6:8.1f} MN/m")
    print(f"  at the hinge, per surface{hinge_stiffness / 1000:8.0f} kNm/rad")
    print()

    print()
    print("DIMENSIONING")

    if force_high > force_low:
        force_case = "high speed"
    else:
        force_case = "low speed"

    if DEFLECTION_HIGH_SPEED > DEFLECTION_LOW_SPEED:
        stroke_case = "high speed"
    else:
        stroke_case = "low speed"

    if speed_high > speed_low:
        speed_case = "high speed"
    else:
        speed_case = "low speed"

    print(f"  rod force   {design_force / 1000:7.1f} kN    from {force_case}")
    print(f"  stroke      {STROKE * 1000:7.1f} mm    from {stroke_case}")
    print(f"  rod speed   {design_speed * 1000:7.1f} mm/s  from {speed_case}")
    print()
    print(
        f"  cylinder    {bore * 1000:.1f} mm bore, {bore * ROD_RATIO * 1000:.1f} mm rod"
    )


if __name__ == "__main__":
    main()
