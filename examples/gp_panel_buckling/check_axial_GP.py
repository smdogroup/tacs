"""
Check axial closed-form modes in GPBladeStiffenedShellConstitutive class are reasonable
@Author Sean Engelstad
@Date 05/16/2024

Use the repo https://github.com/smdogroup/ml_buckling
"""

import numpy as np
import matplotlib.pyplot as plt
import niceplots
from tacs import TACS, constitutive
import ml_buckling as mlb

DEG2RAD = np.pi / 180.0

dtype = TACS.dtype

# build the axial GP object (which is the main ML object we are testing for this example)
# however it is used inside of the constitutive object so we need to build that too
axialGP = constitutive.AxialGP.from_csv(csv_file=mlb.axialGP_csv)
axialGP.setKS(10)
rho0 = 0.3121
xi = 0.9487
gamma = 5.868
zeta = 0.0035

Xtest = np.zeros((4,))
Xtest[0] = np.log(1.0 + xi)
Xtest[1] = np.log(rho0)
Xtest[2] = np.log(1.0 + gamma)
Xtest[3] = np.log(1.0 + 1000.0 * zeta)
pred = axialGP.predict_mean_test_data(Xtest)
print(f"\n", flush=True)
print(f"Xtest = {Xtest}")
print(f"pred = {pred}")
