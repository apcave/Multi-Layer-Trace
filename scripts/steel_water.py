import numpy as np

from modules.composite import Composite

def run_water_steel_water_simulation():
    # Define the properties of the composite
    thickness = [0.0, 1e-3, 0.0]  # Thickness in meters
    density = [1000.0, 7850.0, 1000.0]  # Density in kg/m^3
    cp = [1480.0, 5960.0, 1480.0]  # P-wave speed in m/s
    att_p = [0.01, 0.02, 0.01]  # P-wave attenuation
    cs = [0.0, 3235.0, 0.0]  # S-wave speed in m/s
    att_s = [0.0, 0.05, 0.0]  # S-wave attenuation

    # Create a Composite object
    composite = Composite()
    composite.make_composite(thickness, density, cp, cs, att_p, att_s)

    # Set frequency and angle for the simulation
    frequency = np.linspace(100, 20e3, 100)  # Frequency in Hz
    angle = [30]  # Angle in degrees
    is_compression = True

    composite.set_frequency(frequency)
    composite.set_angle(angle)
    composite.set_is_shear(is_compression)

    # Run the simulation
    composite.run_simulation()

    # Print the results
    if len(composite.fl_rp) <= 10:    
        print("Levesque Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.l_rp)
        print("Transmission coefficient (P-wave):", composite.l_tp)
        #print("Reflection coefficient (S-wave)  :", composite.l_rs)
        #print("Transmission coefficient (S-wave):", composite.l_ts)
        
        print("Cave Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.c_rp)
        print("Transmission coefficient (P-wave):", composite.c_tp)
        #print("Reflection coefficient (S-wave)  :", composite.c_rs)
        #print("Transmission coefficient (S-wave):", composite.c_ts)
        
    composite.plot_results()

def run_water_steel_simulation():
    # Define the properties of the composite
    thickness = [0.0, 0.0]  # Thickness in meters
    density = [1000.0, 7850.0]  # Density in kg/m^3
    cp = [1480.0, 5960.0]  # P-wave speed in m/s
    att_p = [0.01, 0.02]  # P-wave attenuation
    #att_p = [0.0, 0.0]
    #cs = [369.0, 3235.0]  # S-wave speed in m/s
    cs = [369.0, 0.0]  # S-wave speed in m/s
    att_s = [0.1, 0.0]
    att_s = [0.1, 0.05]  # S-wave attenuation
    #att_s = [0.0, 0.0]

    # Create a Composite object
    composite = Composite()
    composite.make_composite(thickness, density, cp, cs, att_p, att_s)

    # Set frequency and angle for the simulation
    frequency = np.linspace(100, 10e3, 100)  # Frequency in Hz
    frequency = [5.0e3]  # Single frequency for the test
    angle = [30]  # Angle in degrees
    is_compression = False

    composite.set_frequency(frequency)
    composite.set_angle(angle)
    composite.set_is_shear(is_compression)

    # Run the simulation
    composite.run_simulation()

    # Print the results
    if len(composite.fl_rp) <= 10:    
        print("Levesque Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.l_rp)
        print("Transmission coefficient (P-wave):", composite.l_tp)
        print("Reflection coefficient (S-wave)  :", composite.l_rs)
        print("Transmission coefficient (S-wave):", composite.l_ts)
        
        print("Cave Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.c_rp)
        print("Transmission coefficient (P-wave):", composite.c_tp)
        print("Reflection coefficient (S-wave)  :", composite.c_rs)
        print("Transmission coefficient (S-wave):", composite.c_ts)
    else:
        composite.plot_results()

def run_steel_water_simulation():
    # Define the properties of the composite
    thickness = [0.0, 0.0]  # Thickness in meters
    density = [7850.0, 1000.0]  # Density in kg/m^3
    cp = [ 5960.0, 1480.0]  # P-wave speed in m/s
    att_p = [0.02, 0.01]  # P-wave attenuation
    #att_p = [0.0, 0.0]
    cs = [3235.0, 0.0]  # S-wave speed in m/s
    att_s = [0.05, 0.0]  # S-wave attenuation
    #att_s = [0.0, 0.0]

    # Create a Composite object
    composite = Composite()
    composite.make_composite(thickness, density, cp, cs, att_p, att_s)

    # Set frequency and angle for the simulation
    frequency = np.linspace(100, 10e3, 100)  # Frequency in Hz
    frequency = [50.0e3]  # Single frequency for the test
    angle = [0]  # Angle in degrees
    is_compression = True

    composite.set_frequency(frequency)
    composite.set_angle(angle)
    composite.set_is_shear(is_compression)

    # Run the simulation
    composite.run_simulation()

    # Print the results
    if len(composite.fl_rp) <= 10:    
        print("Levesque Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.l_rp)
        print("Transmission coefficient (P-wave):", composite.l_tp)
        print("Reflection coefficient (S-wave)  :", composite.l_rs)
        print("Transmission coefficient (S-wave):", composite.l_ts)
        
        print("Cave Acoustic Simulation Results:")
        print("Reflection coefficient (P-wave)  :", composite.c_rp)
        print("Transmission coefficient (P-wave):", composite.c_tp)
        print("Reflection coefficient (S-wave)  :", composite.c_rs)
        print("Transmission coefficient (S-wave):", composite.c_ts)
    else:
        composite.plot_results()

#run_steel_water_simulation()
run_water_steel_simulation()


#run_water_steel_water_simulation()