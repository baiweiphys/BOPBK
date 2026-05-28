function [vdb,vdc] = zeroCurrentCond(v0b,index_c,index_b,par)

% c2 = (2.99792458e8)^2; % speed of ligth c^2
Nb = par.data(index_b,3);
Nc = par.data(index_c,3);


qb = par.data(index_b,1);
qc = par.data(index_c,1);

% Ueb = 1.5e4*1e3; % unit: 1/c

% drift velocity for electron beam (eb)
vdb = v0b; % normalized by 1/c

% zero current condition: Nec*Vec + Neb*Veb = 0
vdc = -Nb*vdb*qb/Nc/qc; % drift velocity for electron core (ec)

disp(['vdb/c = ', num2str(vdb,'%.6e')]);
disp(['vdc/c = ', num2str(vdc,'%.6e')]);
