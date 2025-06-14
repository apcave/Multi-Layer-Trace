import sys, os
sys.path.append(os.path.abspath('./build'))

import composite_wave_response

comp = composite_wave_response.Composite()


# void makeComposite( std::vector<float>& thickness,  std::vector<float>& density,
#                     std::vector<float>& cp,  std::vector<float>& cs,
#                     std::vector<float>& att_p,  std::vector<float>& att_s);

thinkness = [0.0, 5e-3, 0.0]

density = [1000.0, 7850.0, 1000.0]

cp = [1480.0, 5960.0, 1480.0]
att_p = [0.01, 0.02, 0.01]

cs = [ 0.0, 3235.0, 0.0 ]
att_s = [ 0.0, 0.05, 0.0 ]

comp.makeComposite(thinkness, density, cp, cs, att_p, att_s)

