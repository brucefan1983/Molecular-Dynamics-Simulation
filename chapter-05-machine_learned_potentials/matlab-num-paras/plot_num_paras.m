clear; close all;

N_cnk_one = 5 * 13 + 5 * 13;
N_neu = 30;
N_des = 30;
N_ann_one = (N_des + 2) * N_neu 

N_typ = 1:100;
N_par = N_typ * N_ann_one + 1 + N_typ.^2 * N_cnk_one;

figure;
loglog(N_typ, N_par, '-','linewidth',2);
xlabel('$N_{\rm typ}$', 'interpreter', 'latex')
ylabel('$N_{\rm par}$', 'interpreter', 'latex')
set(gca, 'fontsize', 12, 'linewidth', 1)

