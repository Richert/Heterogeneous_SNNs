from pycobi import ODESystem
import sys
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

"""
Bifurcation analysis of the Izhikevich mean-field model with parameters representing an inhibitory fast spiking 
neuron.

To run this code, you need Python >= 3.7 with PyCoBi >= 0.8.4 (https://github.com/pyrates-neuroscience/PyCoBi) 
installed.
"""

# preparations
##############

path = sys.argv[-1]
auto_dir = path if type(path) is str and ".py" not in path else "~/PycharmProjects/auto-07p"

# config
n_dim = 4
n_params = 20
a = ODESystem("fs", working_dir="config", auto_dir=auto_dir, init_cont=False)

# initial continuation in time (to converge to fixed point)
t_sols, t_cont = a.run(c='ivp', name='t', DS=1e-4, DSMIN=1e-10, EPSL=1e-06, NPR=1000, NDIM=n_dim, NPAR=n_params,
                       EPSU=1e-06, EPSS=1e-05, DSMAX=0.1, NMX=50000, UZR={14: 500.0}, STOP={'UZ1'})

########################
# bifurcation analysis #
########################

# main continuation
###################

# continuation in Delta
vals = [0.5]
c1_sols, c1_cont = a.run(starting_point='UZ1', c='qif', ICP=6, NPAR=n_params, NDIM=n_dim, name='D:1',
                         origin=t_cont, NMX=8000, DSMAX=0.01, UZR={6: vals}, STOP=[f'UZ{len(vals)}'], NPR=100,
                         RL1=10.0, RL0=0.0, bidirectional=True)

for i, v in enumerate(vals):

    # continuation in the input strength
    a.run(starting_point=f'UZ{i+1}', c='qif', ICP=16, NPAR=n_params, NDIM=n_dim, name=f'I:{i+1}', origin=c1_cont,
          NMX=8000, DSMAX=0.1, NPR=20, RL1=200.0)

    # continuation of limit cycle
    a.run(starting_point='HB1', c='qif2b', ICP=16, NPAR=n_params, NDIM=n_dim, name=f'I:{i+1}:lc',
          origin=f'I:{i+1}', NMX=10000, DSMAX=0.1, NPR=20, RL1=200.0, STOP=['BP1', 'LP3'])

# 2D continuation follow-up I
target = 0
a.run(starting_point='HB1', c='qif2', ICP=[6, 16], name='D/I:hb1', origin=f'I:{target+1}', NMX=4000, DSMAX=0.2,
      NPR=10, RL1=10.0, RL0=0.001, bidirectional=True)

############
# plotting #
############

fig, ax = plt.subplots(ncols=1, figsize=(6, 6))
a.plot_continuation('PAR(16)', 'PAR(6)', cont=f'D/I:hb1', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta$ (mv)')
ax.set_xlabel(r'$I$ (pA)')
ax.set_title(r'Bifurcation Diagram for Fast Spiking Neurons')
plt.tight_layout()
plt.show()
