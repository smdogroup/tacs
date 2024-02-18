"""
Sean Engelstad
GT SMDO Lab, Feb 2024

Define lambda = load factor to buckling failure
lambda_x = Nxcrit/Nx
lambda_xy = Nxycrit/Nxy
lambda_y = Nycrit/Ny
mu = 1/lambda for all above so that lambda >= 1 equiv to mu <= 1 buckling constraint

Goal is to relate single loading failure criterion
mu_x <= 1; mu_xy <= 1, mu_y <= 1
to combined loading failure criterion
mu = mu(mu_x, mu_y, mu_xy) <= 1
"""

# for a bunch of different designs:
# compute single mode failure criterion by loading only with exx, eyy, or exy
# note that our exx only case actually produces eyy load so that sxx neq 0 is only one
# obtain mu_x, mu_y, mu_xy then (and set mu_i = 0 if that loading style has near zero load)

# compute combined loading failure mode mu
# load mu_x, mu_y, mu_xy, mu and make the failure curve plots
