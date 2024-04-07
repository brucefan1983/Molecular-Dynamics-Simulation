from pylab import *

# Define plotting specifications
aw = 2
fs = 16
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)



# Load in computed and experimental radial distribution functions
aimd_rdf = np.loadtxt("aimd_rdf.out", skiprows=1)
pimd_rdf = np.loadtxt("pimd_rdf.out", skiprows=1)

rdf_exp_goo_T_295_pt1K = np.loadtxt('experimental_rdf/exp-goo-T295.1K.dat', skiprows=6)
rdf_exp_goh_T_300K = np.loadtxt('experimental_rdf/exp-goh-T300K.dat', skiprows=2)
rdf_exp_ghh_T_300K = np.loadtxt('experimental_rdf/exp-ghh-T300K.dat', skiprows=2)

figure()
subplot(3, 1, 1)
set_fig_properties([gca()])
plot(rdf_exp_goo_T_295_pt1K[:, 0], rdf_exp_goo_T_295_pt1K[:, 1], lw=3, ls='--', c='k', label='Experiment')
plot(aimd_rdf[:,0], aimd_rdf[:,3], 'r', lw=2,label='Classical MD')
plot(pimd_rdf[:, 0], pimd_rdf[:, 3], 'b', lw=2, label='Path Integral MD')
ylabel(r'$g_{OO}(r)$', fontsize=22)
xlim([0, 6])
ylim([-0.1, 3])
gca().set_yticks(np.linspace(0, 3, 4))
gca().set_xticklabels([])
xlim([0, 6])
legend(loc=2, fontsize=10, bbox_to_anchor=(0.04, 0.82, 1., .102))

subplot(3, 1,2)
set_fig_properties([gca()])
plot(rdf_exp_goh_T_300K[:, 0], rdf_exp_goh_T_300K[:, 1], lw=3, ls='--', c='k')
plot(aimd_rdf[:,0], aimd_rdf[:,4], 'r', lw=2)
plot(pimd_rdf[:, 0], pimd_rdf[:, 4], 'b', lw=2)
ylabel(r'$g_{OH}(r)$', fontsize=22)
xlim([0, 6])
gca().set_yticks(np.linspace(0, 3, 4))
ylim([-0.1, 2.6])
gca().set_xticks(np.linspace(0, 6, 7))
gca().set_xticklabels([])

subplot(3, 1, 3)
set_fig_properties([gca()])
plot(rdf_exp_ghh_T_300K[:, 0], rdf_exp_ghh_T_300K[:, 1], lw=3, ls='--',c='k')
plot(aimd_rdf[:,0], aimd_rdf[:,2], 'r', lw=2)
plot(pimd_rdf[:, 0], pimd_rdf[:, 2], 'b', lw=2)
ylabel(r'$g_{HH}(r)$', fontsize=22)
xlabel(r'r ($\AA$)', fontsize=22)
xlim([0, 6])
gca().set_yticks(np.linspace(0, 3, 4))
ylim([-0.1, 3.7])
gca().set_xticks(np.linspace(0, 6, 7))
subplots_adjust(hspace=0.1, wspace=0)
savefig('NEP_aimd_pimd_water_rdf.png', dpi=300, bbox_inches='tight', pad_inches=0.1)
show()
