import math
import numpy as np

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
            for ang in self.angle:
                print(f"Running simulation for frequency: {freq} Hz, angle: {ang} degrees, shear: {is_shear}")
                self.cave.properateWave(1.0, 0.0, float(math.radians(ang)),
                                        int(is_shear), float(freq))
                self.fc_rp.append(self.cave.get_rp())
                self.fc_tp.append(self.cave.get_tp())
                self.fc_rs.append(self.cave.get_rs())
                self.fc_ts.append(self.cave.get_ts())
        
    def plot_results(self):
        import matplotlib.pyplot as plt
        
        # Prepare data
        freq_or_angle = self.frequency if len(self.frequency) > 1 else self.angle
        x_label = 'Frequency (Hz)' if len(self.frequency) > 1 else 'Angle (degrees)'

        # Your theory (fc_*)
        rp = np.array(self.fc_rp)
        rs = np.array(self.fc_rs)
        tp = np.array(self.fc_tp)
        ts = np.array(self.fc_ts)
        abs_rp = np.abs(rp)
        abs_rs = np.abs(rs)
        abs_tp = np.abs(tp)
        abs_ts = np.abs(ts)
        phase_rp = np.angle(rp, deg=True)
        phase_rs = np.angle(rs, deg=True)
        phase_tp = np.angle(tp, deg=True)
        phase_ts = np.angle(ts, deg=True)

        # Others (fl_*)
        l_rp = np.array(self.fl_rp)
        l_rs = np.array(self.fl_rs)
        l_tp = np.array(self.fl_tp)
        l_ts = np.array(self.fl_ts)
        abs_l_rp = np.abs(l_rp)
        abs_l_rs = np.abs(l_rs)
        abs_l_tp = np.abs(l_tp)
        abs_l_ts = np.abs(l_ts)
        phase_l_rp = np.angle(l_rp, deg=True)
        phase_l_rs = np.angle(l_rs, deg=True)
        phase_l_tp = np.angle(l_tp, deg=True)
        phase_l_ts = np.angle(l_ts, deg=True)

        fig, axs = plt.subplots(2, 2, figsize=(11.69, 8.27))

        # Top-left: |rp| and |rs|
        axs[0, 0].plot(freq_or_angle, abs_rp, label='|rp| (theory)', color='C0')
        axs[0, 0].plot(freq_or_angle, abs_rs, label='|rs| (theory)', color='C1')
        axs[0, 0].plot(freq_or_angle, abs_l_rp, '--', label='|rp| (other)', color='C0')
        axs[0, 0].plot(freq_or_angle, abs_l_rs, '--', label='|rs| (other)', color='C1')
        axs[0, 0].set_ylabel('Magnitude')
        axs[0, 0].set_xlabel(x_label)
        axs[0, 0].set_title('Reflection Magnitude')
        axs[0, 0].grid(True)
        axs[0, 0].legend()
        axs[0, 0].set_ylim(0, 2)  # Set y-limits for better visibility

        # Bottom-left: |tp| and |ts|
        axs[1, 0].plot(freq_or_angle, abs_tp, label='|tp| (theory)', color='C2')
        axs[1, 0].plot(freq_or_angle, abs_ts, label='|ts| (theory)', color='C3')
        axs[1, 0].plot(freq_or_angle, abs_l_tp, '--', label='|tp| (other)', color='C2')
        axs[1, 0].plot(freq_or_angle, abs_l_ts, '--', label='|ts| (other)', color='C3')
        axs[1, 0].set_ylabel('Magnitude')
        axs[1, 0].set_xlabel(x_label)
        axs[1, 0].set_title('Transmission Magnitude')
        axs[1, 0].grid(True)
        axs[1, 0].legend()
        axs[1, 0].set_ylim(0, 2)  # Set y-limits for better visibility

        # Top-right: phase of rp and rs
        axs[0, 1].plot(freq_or_angle, phase_rp, label='Phase rp (theory)', color='C0')
        axs[0, 1].plot(freq_or_angle, phase_rs, label='Phase rs (theory)', color='C1')
        axs[0, 1].plot(freq_or_angle, phase_l_rp, '--', label='Phase rp (other)', color='C0')
        axs[0, 1].plot(freq_or_angle, phase_l_rs, '--', label='Phase rs (other)', color='C1')
        axs[0, 1].set_ylabel('Phase (degrees)')
        axs[0, 1].set_xlabel(x_label)
        axs[0, 1].set_title('Reflection Phase')
        axs[0, 1].grid(True)
        axs[0, 1].legend()

        # Bottom-right: phase of tp and ts
        axs[1, 1].plot(freq_or_angle, phase_tp, label='Phase tp (theory)', color='C2')
        axs[1, 1].plot(freq_or_angle, phase_ts, label='Phase ts (theory)', color='C3')
        axs[1, 1].plot(freq_or_angle, phase_l_tp, '--', label='Phase tp (other)', color='C2')
        axs[1, 1].plot(freq_or_angle, phase_l_ts, '--', label='Phase ts (other)', color='C3')
        axs[1, 1].set_ylabel('Phase (degrees)')
        axs[1, 1].set_xlabel(x_label)
        axs[1, 1].set_title('Transmission Phase')
        axs[1, 1].grid(True)
        axs[1, 1].legend()

        plt.tight_layout()
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