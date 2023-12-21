from pycobi import ODESystem
import sys
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

"""
Bifurcation analysis of a two-population Izhikevich mean-field model with an excitatory regular spiking neuron population
and an inhibitory fast spiking neuron population.

To run this code, you need Python >= 3.6 with PyCoBi >= 0.8.4 (https://github.com/pyrates-neuroscience/PyCoBi) 
installed.
"""

# preparations
##############

path = sys.argv[-1]
auto_dir = path if type(path) is str and ".py" not in path else "~/PycharmProjects/auto-07p"

# config
n_dim = 8
n_params = 33
a = ODESystem('eic_shadowing', working_dir="config", auto_dir=auto_dir, init_cont=False)

# initial continuation in time (to converge to fixed point)
t_sols, t_cont = a.run(c='ivp', name='t', DS=1e-4, DSMIN=1e-12, EPSL=1e-06, EPSU=1e-06, EPSS=1e-04,
                       DSMAX=0.1, NMX=50000, UZR={14: 2000.0}, STOP={'UZ1'}, NPR=1000, NDIM=n_dim, NPAR=n_params)

########################
# bifurcation analysis #
########################

# Set up the conditions
#######################

# continuation in background input to RS population
c1_sols, c1_cont = a.run(starting_point='UZ1', c='qif', ICP=15, NPAR=n_params, NDIM=n_dim, name='I_rs:1',
                         origin=t_cont, NMX=8000, DSMAX=0.1, UZR={15: [60.0]}, STOP=[f'UZ1'], NPR=100,
                         RL1=100.0)

# continuation in Delta_fs
vals = [0.2, 2.0]
c2_sols, c2_cont = a.run(starting_point='UZ1', c='qif', ICP=24, NPAR=n_params, NDIM=n_dim, name='D_fs:1',
                         origin=c1_cont, NMX=8000, DSMAX=0.1, UZR={24: vals}, STOP=[], NPR=100, RL1=3.0,
                         RL0=0.0, bidirectional=True)

# continuation in d_rs for low Delta_fs
vals = [10.0, 100.0]
c3_sols, c3_cont = a.run(starting_point='UZ1', c='qif', ICP=18, NPAR=n_params, NDIM=n_dim, name='d_rs:1',
                         origin=c2_cont, NMX=8000, DSMAX=0.1, UZR={18: vals}, STOP=[], NPR=100, RL1=150.0,
                         RL0=0.0, bidirectional=True)

# continuation in d_rs for high Delta_fs
vals = [10.0, 100.0]
c4_sols, c4_cont = a.run(starting_point='UZ2', c='qif', ICP=18, NPAR=n_params, NDIM=n_dim, name='d_rs:2',
                         origin=c2_cont, NMX=8000, DSMAX=0.1, UZR={18: vals}, STOP=[], NPR=100, RL1=150.0,
                         RL0=0.0, bidirectional=True)

# 2D bifurcation analysis in Delta_rs and I_fs
##############################################

# continuation in I_fs for Delta_fs = 0.2 and d_rs = 10.0
r1_sols, r1_cont = a.run(starting_point='UZ1', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:1',
                         origin=c3_cont, NMX=8000, DSMAX=0.1, UZR={26: [45.0]}, STOP=[], NPR=100, RL1=150.0)
a.run(starting_point='HB1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:1:hb1', origin=r1_cont,
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2', 'BP1'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)
a.run(starting_point='HB3', c='qif2', ICP=[6, 30], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:1:hb2', origin=r1_cont,
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)

# continuation in I_fs for Delta_fs = 0.2 and d_rs = 100.0
r2_sols, r2_cont = a.run(starting_point='UZ2', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:2',
                         origin=c3_cont, NMX=8000, DSMAX=0.1, UZR={}, STOP=[], NPR=100, RL1=150.0)
a.run(starting_point='HB1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:2:hb1', origin=r2_cont,
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)
a.run(starting_point='HB2', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:2:hb2', origin=r2_cont,
      NMX=8000, DSMAX=0.1, UZR={6: [0.122]}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)
a.run(starting_point='UZ1', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:2:2', origin='D_rs/I_fs:2:hb2',
      NMX=8000, DSMAX=0.1, UZR={}, STOP=[], NPR=100, RL0=0.0, bidirectional=True)
a.run(starting_point='LP1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:2:lp1', origin='I_fs:2:2',
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)

# continuation in I_fs for Delta_fs = 2.0 and d_rs = 10.0
r3_sols, r3_cont = a.run(starting_point='UZ1', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:3',
                         origin=c4_cont, NMX=8000, DSMAX=0.1, UZR={}, STOP=[], NPR=100, RL1=150.0)
a.run(starting_point='LP1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:3:lp1', origin=r3_cont,
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)
a.run(starting_point='LP2', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:3:lp2', origin=r3_cont,
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)

# continuation in I_fs for Delta_fs = 2.0 and d_rs = 100.0
r4_sols, r4_cont = a.run(starting_point='UZ2', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:4',
                         origin=c4_cont, NMX=8000, DSMAX=0.1, UZR={}, STOP=[], NPR=100, RL1=200.0)
a.run(starting_point='HB1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:4:hb1', origin=r4_cont,
      NMX=8000, DSMAX=0.1, UZR={26: [27.75]}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)
a.run(starting_point='UZ1', c='qif', ICP=6, NPAR=n_params, NDIM=n_dim, name='D_rs:4', origin='D_rs/I_fs:4:hb1',
      NMX=8000, DSMAX=0.1, UZR={6: [0.2]}, STOP=["UZ1"], NPR=100, RL0=0.0, DS="-")
a.run(starting_point='UZ1', c='qif', ICP=26, NPAR=n_params, NDIM=n_dim, name='I_fs:4:2', origin='D_rs:4',
      NMX=8000, DSMAX=0.1, UZR={}, STOP=[], NPR=100, RL1=100.0, RL0=0.0, bidirectional=True)
a.run(starting_point='LP1', c='qif2', ICP=[6, 26], NPAR=n_params, NDIM=n_dim, name='D_rs/I_fs:4:lp1', origin='I_fs:4:2',
      NMX=8000, DSMAX=0.1, UZR={}, STOP=['CP2'], NPR=10, RL1=10.0, RL0=0.0, bidirectional=True, EPSS=1e-6)

############
# plotting #
############

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

# Delta_fs = 1.0, d_rs = 10.0
ax = axes[0, 0]
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:3:lp1', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:3:lp2', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta_{rs}$ (mv)')
ax.set_xlabel(r'$I_{fs}$ (pA)')
ax.set_title(r'(A) $\Delta_{fs} = 2.0$ mV, $\kappa_{rs} = 10.0$')
ax.set_ylim([0.0, 2.5])
ax.set_xlim([70.0, 10.0])

# Delta_fs = 1.0, d_rs = 100.0
ax = axes[0, 1]
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:4:hb1', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:4:lp1', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta_{rs}$ (mv)')
ax.set_xlabel(r'$I_{fs}$ (pA)')
ax.set_title(r'(B) $\Delta_{fs} = 2.0$ mV, $\kappa_{rs} = 100.0$')
ax.set_ylim([0.0, 2.5])
ax.set_xlim([70.0, 10.0])

# Delta_fs = 0.1, d_rs = 10.0
ax = axes[1, 0]
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:1:hb1', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
a.plot_continuation('PAR(30)', 'PAR(6)', cont=f'D_rs/I_fs:1:hb2', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta_{rs}$ (mv)')
ax.set_xlabel(r'$I_{fs}$ (pA)')
ax.set_title(r'(C) $\Delta_{fs} = 0.2$ mV, $\kappa_{rs} = 10.0$')
ax.set_ylim([0.0, 2.5])
ax.set_xlim([70.0, 10.0])

# EIC: Delta_fs = 0.1, d_rs = 100.0
ax = axes[1, 1]
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:2:hb1', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:2:hb2', ax=ax, line_color_stable='#148F77',
                    line_color_unstable='#148F77', line_style_unstable='solid')
a.plot_continuation('PAR(26)', 'PAR(6)', cont=f'D_rs/I_fs:2:lp1', ax=ax, line_color_stable='#5D6D7E',
                    line_color_unstable='#5D6D7E', line_style_unstable='solid')
ax.set_ylabel(r'$\Delta_{rs}$ (mv)')
ax.set_xlabel(r'$I_{fs}$ (pA)')
ax.set_title(r'(D) $\Delta_{fs} = 0.2$ mV, $\kappa_{rs} = 100.0$')
ax.set_ylim([0.0, 2.5])
ax.set_xlim([70.0, 10.0])

plt.suptitle("Bifurcation Diagrams for a Two-Population Model of Regular Spiking (RS) and Fast Spiking (FS) Neurons")
plt.tight_layout()
plt.show()
