from pycobi import ODESystem
import sys
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

"""
Bifurcation analysis of the Izhikevich mean-field model with parameters representing an excitatory regular spiking 
neuron.

To run this code, you need Python >= 3.7 with PyCoBi >= 0.8.4 (https://github.com/pyrates-neuroscience/PyCoBi) 
installed.
"""

# preparations
##############

# directory of the auto-07p installation
path = sys.argv[-1]
auto_dir = path if type(path) is str and ".py" not in path else "~/PycharmProjects/auto-07p"

# config
n_dim = 4
n_params = 20
a = ODESystem("rs", working_dir="config", auto_dir=auto_dir, init_cont=False)

# initial continuation in time to converge to fixed point
t_sols, t_cont = a.run(c='ivp', name='t', DS=1e-4, DSMIN=1e-10, EPSL=1e-06, NPR=1000, NPAR=n_params, NDIM=n_dim,
                       EPSU=1e-06, EPSS=1e-05, DSMAX=0.1, NMX=50000, UZR={14: 500.0}, STOP={'UZ1'})

########################
# bifurcation analysis #
########################

# prepare state
###############

# continuation in synaptic strength
c1_sols, c1_cont = a.run(starting_point='UZ1', c='qif', ICP=4, NPAR=n_params, NDIM=n_dim, name='g:1',
                         origin=t_cont, NMX=8000, DSMAX=0.05, UZR={4: [15.0]}, STOP=[f'UZ1'], NPR=100,
                         RL1=1000.0, RL0=0.0)

# continuation in Delta
vals1 = [0.1, 1.0]
c2_sols, c2_cont = a.run(starting_point='UZ1', c='qif', ICP=5, NPAR=n_params, NDIM=n_dim, name='D:1',
                         origin=c1_cont, NMX=8000, DSMAX=0.05, UZR={5: vals1}, STOP=[f'UZ{len(vals1)}'], NPR=100,
                         RL1=3.0, RL0=0.0, bidirectional=True)

# continuations in d
vals = [10.0, 100.0]
c3_sols, c3_cont = a.run(starting_point='UZ1', c='qif', ICP=16, NPAR=n_params, NDIM=n_dim, name='d:1',
                         origin=c2_cont, NMX=8000, DSMAX=0.05, UZR={16: vals}, STOP=[f'UZ{len(vals)}'], NPR=100,
                         RL1=150.0, RL0=0.0, bidirectional=True)
c4_sols, c4_cont = a.run(starting_point='UZ2', c='qif', ICP=16, NPAR=n_params, NDIM=n_dim, name='d:2',
                         origin=c2_cont, NMX=8000, DSMAX=0.05, UZR={16: vals}, STOP=[f'UZ{len(vals)}'], NPR=100,
                         RL1=150.0, RL0=0.0, bidirectional=True)

# main continuations
####################

# continuations in I for d = 10.0
r1_sols, r1_cont = a.run(starting_point='UZ1', c='qif', ICP=8, NPAR=n_params, NDIM=n_dim, name='I_ext:1',
                         origin=c3_cont, NMX=8000, DSMAX=0.05, UZR={}, STOP=[], NPR=100,
                         RL1=150.0, RL0=0.0, bidirectional=True)
a.run(starting_point='LP1', c='qif2', ICP=[5, 8], name='D/I:lp1', origin=r1_cont, NMX=8000, DSMAX=0.05,
      NPR=20, RL1=5.0, RL0=0.0, bidirectional=True)
a.run(starting_point='LP2', c='qif2', ICP=[5, 8], name='D/I:lp2', origin=r1_cont, NMX=8000, DSMAX=0.05,
      NPR=20, RL1=5.0, RL0=0.0, bidirectional=True)

# continuations in I for d = 100.0
r2_sols, r2_cont = a.run(starting_point='UZ2', c='qif', ICP=8, NPAR=n_params, NDIM=n_dim, name='I_ext:2',
                         origin=c3_cont, NMX=8000, DSMAX=0.05, UZR={}, STOP=[], NPR=100,
                         RL1=150.0, RL0=0.0, bidirectional=True)
r3_sols, r3_cont = a.run(starting_point='UZ2', c='qif', ICP=8, NPAR=n_params, NDIM=n_dim, name='I_ext:3',
                         origin=c4_cont, NMX=8000, DSMAX=0.05, UZR={}, STOP=[], NPR=100,
                         RL1=150.0, RL0=0.0, bidirectional=True)
a.run(starting_point='LP1', c='qif2', ICP=[5, 8], name='D/I:lp3', origin=r2_cont, NMX=8000, DSMAX=0.05,
      NPR=20, RL1=5.0, RL0=0.0, bidirectional=True)
a.run(starting_point='LP2', c='qif2', ICP=[5, 8], name='D/I:lp4', origin=r2_cont, NMX=8000, DSMAX=0.05,
      NPR=20, RL1=5.0, RL0=0.0, bidirectional=True)
a.run(starting_point='HB1', c='qif2', ICP=[5, 8], name='D/I:hb1', origin=r3_cont, NMX=8000, DSMAX=0.05,
      NPR=20, RL1=5.0, RL0=0.0, bidirectional=True)

############
# plotting #
############

fig, axes = plt.subplots(ncols=2, figsize=(12, 6))

ax = axes[0]
a.plot_continuation('PAR(8)', 'PAR(5)', cont=f'D/I:lp1', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
a.plot_continuation('PAR(8)', 'PAR(5)', cont=f'D/I:lp2', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta$ (mv)')
ax.set_xlabel(r'$I$ (pA)')
ax.set_title(r'(A) $\kappa = 10.0$ pA')

ax = axes[1]
a.plot_continuation('PAR(8)', 'PAR(5)', cont=f'D/I:lp3', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
a.plot_continuation('PAR(8)', 'PAR(5)', cont=f'D/I:lp4', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
a.plot_continuation('PAR(8)', 'PAR(5)', cont=f'D/I:hb1', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta$ (mv)')
ax.set_xlabel(r'$I$ (pA)')
ax.set_title(r'(B) $\kappa = 100.0$ pA')
plt.suptitle("Bifurcation Diagram for Regular Spiking Neurons")
plt.tight_layout()
plt.show()
