%% Mean field theory (MFT) parameters form DyCoNi
% Author: Rasmus Bjørk, Technical University of Denmark, rabj@dtu.dk
%
%--- The scripts shows how to use the MFT model, with input parameters for the material DyCoNi
%
%--- Please cite these works as well as this repository when utilizing the model to generate data.

clearvars

Tc      = 67.1;                                 % Curie temperature
gj      = 1.33;                                 % Landé factor
J       = 7.5;                                  % Total angular momentum
thetaD  = 350;                                  % Debye temperature (K)
Ns      = 1.2594e+24;                           % Number of spins pr unit mass in kg^-1
M       = 1/3*(162.5e-3+58.933e-3+58.693e-3);   % Molar mass [kg/mol]
gamma_e = 0;                                    % Sommerfeld constant (J/kgK^2)

Tarr = linspace(2,100,20);
Barr = [5e-5 linspace(1,5,5)];

x0 = [Tc gj J thetaD Ns M gamma_e];
[Cp_tot, DT, S, DS, mag]  = MFT_model(x0,Tarr,Barr,'SaveTheResult',false,'ShowTheResult',false);

%% ------------------------------------------------------------------
%% -------------------- Visualize the results -----------------------
%% ------------------------------------------------------------------
colorarr = turbo(length(Cp_tot(:,1)));

figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); hold on; grid on;
for i = 1:length(Cp_tot(:,1))
    plot(Tarr,Cp_tot(i,:),'Marker','.','Markersize',10,'Linestyle','none','Color',colorarr(i,:));
end
xlim([min(Tarr) max(Tarr)]);
xlabel('Temperature [K]')
ylabel('c_p [J kg^{-1} K^{-1}]')

figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); hold on; grid on;
for i = 1:length(DT(:,1))
    plot(Tarr,DT(i,:),'Marker','.','Markersize',10,'Linestyle','none','Color',colorarr(i,:));
end
xlabel('Temperature [K]')
ylabel('\Delta{}T_{ad} [K]')

figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); hold on; grid on;
for i = 1:length(DS(:,1))
    plot(Tarr,DS(i,:),'Marker','.','Markersize',10,'Linestyle','none','Color',colorarr(i,:));
end
xlabel('Temperature [K]')
ylabel('\Delta{}s_{iso} [J kg^{-1} K^{-1}]')

figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig4 = axes('Parent',figure4,'Layer','top','FontSize',16); hold on; grid on;
for i = 1:length(S(:,1))
    plot(Tarr,S(i,:),'Marker','.','Markersize',10,'Linestyle','none','Color',colorarr(i,:));
end
xlabel('Temperature [K]')
ylabel('S')

figure5= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig5 = axes('Parent',figure5,'Layer','top','FontSize',16); hold on; grid on;
for i = 1:length(mag(:,1))
    plot(Tarr,mag(i,:),'Marker','.','Markersize',10,'Linestyle','none','Color',colorarr(i,:));
end
xlabel('Temperature [K]')
ylabel('\sigma [emu/g or Am^2 kg^{-1}]')
