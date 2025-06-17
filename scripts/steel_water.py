import numpy as np

from modules.composite import Composite

composite = Composite()


def make_water_steel_water():
    # Define the properties of the composite
    thickness = [0.0, 30e-3, 0.0]  # Thickness in meters
    density = [1000.0, 900, 1000.0]  # Density in kg/m^3
    cp = [1480.0, 1260.0, 1480.0]  # P-wave speed in m/s

    att_p = [0.01, 0.2, 0.01]  # P-wave attenuation
    cs = [0.0, 337, 0.0]  # S-wave speed in m/s
    att_s = [0.0, 1, 0.0]  # S-wave attenuation

    composite.make_composite(thickness, density, cp, cs, att_p, att_s)

def make_water_steel():
    # Define the properties of the composite
    thickness = [0.0, 0.0]  # Thickness in meters
    density = [1000.0, 7850.0]  # Density in kg/m^3
    cp = [1480.0, 5960.0]  # P-wave speed in m/s
    att_p = [0.01, 0.02]  # P-wave attenuation
    cs = [0, 3235.0]  # S-wave speed in m/s
    att_s = [0, 0.1]
    composite.make_composite(thickness, density, cp, cs, att_p, att_s)


def make_steel_water():
    # Define the properties of the composite
    thickness = [0.0, 0.0]  # Thickness in meters
    density = [7850.0, 1000.0]  # Density in kg/m^3
    cp = [ 5960.0, 1480.0]  # P-wave speed in m/s
    att_p = [0.02, 0.01]  # P-wave attenuation
    #att_p = [0.0, 0.0]
    cs = [3235.0, 0.0]  # S-wave speed in m/s
    att_s = [0.05, 0.0]  # S-wave attenuation
    #att_s = [0.0, 0.0]

    composite.make_composite(thickness, density, cp, cs, att_p, att_s)


def make_water_fluid():
    # Define the properties of the composite
    thickness = [0.0, 1e-3, 0.0]  # Thickness in meters
    density = [1000.0, 7850.0, 1000.0]  # Density in kg/m^3
    cp = [1480.0, 5960.0, 1480.0]  # P-wave speed in m/s
    att_p = [1, 20, 1]  # P-wave attenuation
    cs = [0.0, 0, 0.0]  # S-wave speed in m/s
    att_s = [0.0, 0, 0.0]  # S-wave attenuation

    #att_p = [0,0,0]
    composite.make_composite(thickness, density, cp, cs, att_p, att_s)

def make_solid_solid_solid():
    # Define the properties of the composite
    thickness = [0.0, 0.0, 0.0]  # Thickness in meters
    density = [7850.0, 1000, 7850.0]  # Density in kg/m^3
    cp = [5960.0, 1480, 5960.0]  # P-wave speed in m/s
    att_p = [1, 3, 1]  # P-wave attenuation
    cs = [3235.0, 413, 3235.0]  # S-wave speed in m/s
    att_s = [1, 7.9, 1]  # S-wave attenuation

    composite.make_composite(thickness, density, cp, cs, att_p, att_s)


# make_water_steel_water()
# make_water_fluid()
# make_solid_solid_solid()
make_water_steel_water()
# Set frequency and angle for the simulation
frequency = np.linspace(100, 20e3, 101)  # Frequency in Hz
frequency = [10.0e3]  # Single frequency for the test
angle = np.linspace(0, 90, 91)  # Angle in degrees
#angle = [15]  # Angle in degrees
is_compression = True

composite.set_frequency(frequency)
composite.set_angle(angle)
composite.set_is_shear(is_compression)

# Run the simulation
composite.run_simulation()

# Print the results
if len(composite.fl_rp) <= 10:    
    print("Levesque Acoustic Simulation Results:")
    print("Reflection coefficient (P-wave)  :", abs(composite.l_rp), np.angle(composite.l_rp))
    print("Transmission coefficient (P-wave):", abs(composite.l_tp), np.angle(composite.l_tp))
    print("Reflection coefficient (S-wave)  :", abs(composite.l_rs), np.angle(composite.l_rs))
    print("Transmission coefficient (S-wave):", abs(composite.l_ts), np.angle(composite.l_ts))
    
    print("Cave Acoustic Simulation Results:")
    print("Reflection coefficient (P-wave)  :", abs(composite.c_rp), np.angle(composite.c_rp))
    print("Transmission coefficient (P-wave):", abs(composite.c_tp), np.angle(composite.c_tp))
    print("Reflection coefficient (S-wave)  :", abs(composite.c_rs), np.angle(composite.c_rs))
    print("Transmission coefficient (S-wave):", abs(composite.c_ts), np.angle(composite.c_ts))
else:
    composite.plot_results()