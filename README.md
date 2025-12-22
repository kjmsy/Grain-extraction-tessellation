
Tessellation-based grain boundary determination for X-ray orientation microscopies

Jaemyung Kim, Yujiro Hayashi, Sung Soo Ha, Makina Yabashi,
Computer Physics Communications (2026). 320, 109988
https://doi.org/10.1016/j.cpc.2025.109988.
(https://www.sciencedirect.com/science/article/pii/S0010465525004898)

Figures in the paper can be reproduced from the Python scripts.


The code is developed under Python 3.13.2 together with Spyder.

The input variables in the files are 
data file name, misorientation threshold, IPF axis, and crystalline phase.


1. The header of the input data file (example_1.txt) is 

x y u_11 u_12 u_13 u_21 u_22 u_23 u_31 u_32 u_33

x, and y are position.

u_11 u_12 u_13 u_21 u_22 u_23 u_31 u_32 u_33 --> Orientation matrix U

U = | u_11 u_12 u_13 |

    | u_21 u_22 u_23 |

    | u_31 u_32 u_33 |

The transpose of U is g, which is usually used in the EBSD (Electron back scatter diffraction) measurement.
    

2. Misorientation threshold
Let's consider 2 neighboring pixels.
If their misorientations are below the threshold, we can consider them as a connected grain.
If not, they are separated grain.
This value should be larger than 0.
The value is the degree unit in the script.


3. IPF axis
There are 3 axes in the inverse pole figure (IPF).
Metallurgists describe them as RD (rolling direction), TD (transverse direction), and ND (normal direction).
They are IPF-X, IPF-Y, and IPF-Z, respectively.
They should be "ipf_x", "ipf_y", and "ipf_z".


4. phase
The IPF is dependent on the crystal symmetry.
Therefore, the crystalline phase should be defined for the IPF plot.
This version supports cubic, hexagonal, and tetragonal.
To set it to cubic, set "phase = 0".
For hexagonal and tetragonal, set it to "1", and "2", respectively.

