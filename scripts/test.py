import math
import numpy as np
import sys, os
sys.path.append(os.path.abspath('./build'))

from composite_wave_response import makeComposite, properateWave


# void makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
#                     std::vector<float>& cp,  std::vector<float>& cs,
#                     std::vector<float>& att_p,  std::vector<float>& att_s);


thickness = [0.0, 5e-3, 0.0]

density = [1000.0, 7850.0, 1000.0]

cp =      [1480.0, 5960.0, 1480.0]
att_p = [0.01, 0.02, 0.01]

cs =    [ 0.0, 3235.0, 0.0 ]
att_s = [ 0.0, 0.05, 0.0 ]

angle = 0.0
frequency = 1.0e3
isShear = False

makeComposite(thickness, density, cp, cs, att_p, att_s)

properateWave(1.0, 0.0, float(math.radians(angle)),
              int(isShear), float(frequency))