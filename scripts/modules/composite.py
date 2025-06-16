import math

from modules import ml_cave
from modules import ml_levesque

class Composite:
    def __init__(self):
        self.levesque = ml_levesque
        self.cave = ml_cave
        self.fl_rp = []
        self.fl_tp = []
        self.fl_rs = []
        self.fl_ts = []

    def make_composite(self, thickness, density, cp, cs, att_p, att_s):
        number_mediums = len(thickness)
        if number_mediums < 2:
            raise ValueError("At least two mediums are required to create a composite.")
        if len(density) != number_mediums or len(cp) != number_mediums or len(cs) != number_mediums or len(att_p) != number_mediums or len(att_s) != number_mediums:
            raise ValueError("All input lists must have the same length as the number of mediums.")
        self.cave.makeComposite(thickness, density, cp, cs, att_p, att_s)
        self.thickness = thickness
        self.density = density
        self.cp = cp
        self.cs = cs
        self.att_p = att_p
        self.att_s = att_s
        self.LongM = [0.0] * number_mediums
        self.mu = [0.0] * number_mediums
        
    def set_frequency(self, frequency):
        self.frequency = frequency
        
    def set_angle(self, angle):
        self.angle = angle
        
    def set_is_shear(self, is_compression):
        self.is_compression = is_compression
        
    def run_simulation(self):
        if not hasattr(self, 'frequency') or not hasattr(self, 'angle') or not hasattr(self, 'is_compression'):
            raise ValueError("Frequency, angle, and is_compression must be set before running the simulation.")
        
        if len(self.frequency) > 1 and len(self.angle) > 1:
            raise ValueError("Only one frequency and one angle can be set for the simulation.")
             
        self.fl_rp, self.fl_tp, self.fl_rs, self.fl_ts = self.levesque.run_acoustic_simulation(
            self.is_compression, self.angle, self.frequency, self.density,
            self.thickness, self.cp, self.att_p,
            self.LongM, self.cs, self.att_s, self.mu
        )
        
        self.fc_rp = []
        self.fc_tp = []
        self.fc_rs = []
        self.fc_ts = []
        
        is_shear = 0 if self.is_compression else 1
        for freq in self.frequency:
            for angle in self.angle:
                self.cave.properateWave(1.0, 0.0, float(math.radians(angle)),
                                        int(is_shear), float(freq))
                self.fc_rp.append(self.cave.get_rp())
                self.fc_tp.append(self.cave.get_tp())
                self.fc_rs.append(self.cave.get_rs())
                self.fc_ts.append(self.cave.get_ts())
        
    def plot_results(self):
        import matplotlib.pyplot as plt
        
        if len(self.frequency) > 1:
            plt.figure(figsize=(10, 6))
            plt.plot(self.frequency, [abs(r) for r in self.fc_rp], label='Cave Rp')
            plt.plot(self.frequency, [abs(t) for t in self.fc_tp], label='Cave Tp')
            plt.plot(self.frequency, [abs(r) for r in self.fl_rp], label='Levesque Rp')
            plt.plot(self.frequency, [abs(t) for t in self.fl_tp], label='Levesque Tp')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Magnitude')
            plt.title('Levesque Acoustic Simulation Results')
            plt.legend()
            plt.grid(True)
            plt.show()
        
        if len(self.angle) > 1:
            plt.figure(figsize=(10, 6))
            plt.plot(self.angle, [abs(r) for r in self.fc_rp], label='Reflection Coefficient (P-wave)')
            plt.plot(self.angle, [abs(t) for t in self.fc_tp], label='Transmission Coefficient (P-wave)')
            plt.xlabel('Angle (degrees)')
            plt.ylabel('Magnitude')
            plt.title('Cave Acoustic Simulation Results')
            plt.legend()
            plt.grid(True)
            plt.show()
        
    @property        
    def l_rp(self):
        if len(self.fl_rp) == 1:
            return self.fl_rp[0]
        return self.fl_rp
    
    @property
    def l_tp(self):
        if len(self.fl_tp) == 1:
            return self.fl_tp[0]
        return self.fl_tp
    
    @property
    def l_rs(self):
        if len(self.fl_rs) == 1:
            return self.fl_rs[0]
        return self.fl_rs
    
    @property
    def l_ts(self):
        if len(self.fl_ts) == 1:
            return self.fl_ts[0]
        return self.fl_ts
    
    @property
    def c_rp(self):
        if len(self.fc_rp) == 1:
            return self.fc_rp[0]
        return self.fc_rp
    
    @property
    def c_tp(self):
        if len(self.fc_tp) == 1:
            return self.fc_tp[0]
        return self.fc_tp
    
    @property
    def c_rs(self):
        if len(self.fc_rs) == 1:
            return self.fc_rs[0]
        return self.fc_rs
    
    @property
    def c_ts(self):
        if len(self.fc_ts) == 1:
            return self.fc_ts[0]
        return self.fc_ts