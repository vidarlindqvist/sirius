# Sirius — rudder actuator sizing

Sizes the hydraulic actuator for the split rudder from the SAAB duty points,
through the linkage, to a pump specification.

    python sizing.py

## What it produces

    FLOW AND POWER    the pump specification
    PUSHROD LOAD      the largest force in the mechanism
    STIFFNESS         for the bandwidth and flutter assessment
    DIMENSIONING      which flight case sizes what


## Inputs

Everything you can change is at the top of `sizing.py`, in four groups.

### Duty — from SAAB

| Constant | Meaning |
|---|---|
| `MOMENT_HIGH_SPEED`, `RATE_*`, `DEFLECTION_*` | Mach 1.6, q 47 kPa |
| `MOMENT_LOW_SPEED`, `RATE_*`, `DEFLECTION_*` | Mach 0.2, q 3 kPa |

Hinge moments come from CFD. Rates and deflections are the requirement.

### Architecture

| Constant | Meaning |
|---|---|
| `MECHANISMS` | actuators per surface, sharing the hinge moment |
| `SURFACES` | upper + lower, left + right |
| `LOAD_PATHS` | ball links sharing the load at the horn |

### From the linkage

These five are the only inputs that come from the **mechanism** rather than
from a requirement. Get them by driving the rudder through its travel in CAD,
logging actuator length against rudder angle, and taking the slope.

| Constant | How to get it |
|---|---|
| `LEVER_ARM_HIGH_SPEED` | slope of that curve at the high-speed deflection |
| `LEVER_ARM_LOW_SPEED` | slope at the low-speed deflection |
| `LEVER_ARM_PEAK` | largest slope anywhere in the travel |
| `STROKE` | longest minus shortest actuator length |
| `HORN_LENGTH` | hinge to pushrod pin, straight off the drawing |

The lever arm is a Jacobian — metres of rod travel per radian of rudder.
Everything downstream follows from it.

### Hydraulics

`SUPPLY_PRESSURE`, `ROD_RATIO`, `EFFICIENCY`, `BULK_MODULUS`.


## Reading the output

**Read DIMENSIONING first.** It reports three separate maxima — force, stroke
and speed — and says which flight case won each. They need not be the same
case, and usually are not. Sizing everything from a single "worst case" gets
the stroke badly wrong.

**In FLOW AND POWER**, the cap-side figure is the biggest flow in the circuit.
It comes from the pump when extending and goes to tank when retracting, so
return lines must pass it too.

**In PUSHROD LOAD**, note that the total is divided by `LOAD_PATHS`. Pin
diameter goes as the square root of load, so four joints are half the
diameter, not a quarter — that is how the joint is made to fit inside the
wing thickness.

**STIFFNESS is an upper bound.** Hoses, line volume, back-up structure and
entrained air all reduce it. Nothing increases it.


## A built-in check

Flow is also computed a second way, from a closed form containing no lever
arm and no bore:

    Q = M * rate / (MECHANISMS * EFFICIENCY * PRESSURE * (1 - ROD_RATIO^2))

The two should agree to a few percent. If they do not, one of exactly two
things is wrong: the cylinder is larger than the load requires, or the lever
arm collapses somewhere in the travel.


## Units

SI throughout — metres, pascals, newtons, radians. Conversion to mm, kN and
L/min happens only in the print statements. Do not mix units inside the file.


## What this does NOT do

  - solve the linkage, or check it for interference and singularities
  - any structural check: no wall thickness, buckling, pin, lug or bearing
  - fatigue — everything here is static
  - the bellcrank body and horn root, which need FEA

It takes the lever arms as given and does the sizing arithmetic that follows.


## Provisional

Lever arms and horn length come from a mockup linkage, not from released CAD.
`EFFICIENCY`, `ROD_RATIO` and `BULK_MODULUS` are chosen, not derived. Nothing
has been verified against hardware.
