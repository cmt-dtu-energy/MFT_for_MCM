%clearvars

addpath('..\..\MFT_for_MCM\')

rho = 7901; % Density in kg/m^3
k   = 5;    % Thermal conductivity in W/m-K from Watanabe 2021.

Tc      = 293.6;                 % Curie temperature [K]
gj      = 2;                  % Landé factor
J       = 3.5;                % Total angular momentum
thetaD  = 169;                % Debye temperature (K)
Ns      = 3.83*1e24;          % Number of spins pr unit mass in kg^-1
M       = 0.15725;            % Molar mass [kg/mol]
gamma_e = 6.93*1e-2;           % Sommerfeld constant (J/kgK^2)


x0 = [Tc, gj, J, thetaD, Ns, M, gamma_e];

Tarr = linspace(253,353,101);
Barr = [5e-5 linspace(0.05,3,60)];
[Cp_tot, DT, S, DS, mag]  = MFT_model(x0,Tarr,Barr,'SaveTheResult',false,'ShowTheResult',true);

save('Gd_MFT.mat','Barr','DS','Tarr','S','mag','x0','ci', 'rho', 'k');