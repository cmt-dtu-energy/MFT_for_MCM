function [Cp, DT, S, DS, mag]  = MFT_model(x0,Tarr,Barr,options)
%% Mean field theory (MFT) model for magnetocalorics
% Author: Rasmus Bjørk, Technical University of Denmark, rabj@dtu.dk

%--- This repository contains an implementation in Matlab of the mean field theory (MFT) model, also called the mean field model (MFM), for generating material data for magnetocaloric materials (MCMs).

%--- The program takes as input the Curie temperature, the Landé factor, the total angular momentum, the Debye temperature, the number of spins pr unit mass, the molar mass and the Sommerfeld constant.
%--- For a specified temperature range and magnetic field range, the program outputs the specific heat, the entropy, the magnetization and the adiabatic temperature change of the material.

%--- The theory of the model can be found on p.6 (sec. 2.2.1) in the publication https://orbit.dtu.dk/files/101639258/Designing_a_magnet.pdf 
%--- or in [Smith, A., Bahl, C.R.H., Bjørk, R., Engelbrecht, K., Nielsen, K. K. and Pryds, N., ...
%--- "Materials challenges for high performance magnetocaloric refrigeration devices." ...
%--- Advanced Energy Materials 2, no. 11 (2012): 1288-1318, DOI: 10.1002/aenm.201200167. 
%
%--- Please cite these works as well as this repository when utilizing the model to generate data.

%%
arguments
    x0 (1,7) {mustBeNumeric}                = [298 2 3.5 169 2.7816*1e24 0.15725 0.0109/0.15725]   %--- Material properties for Gd
    Tarr (1,:) {mustBeNumeric}              = 200:1:400;     %--- Temperature values
    Barr  (1,:) {mustBeNumeric}             = [0.01 1 5]     %--- Field values
    options.ShowTheResult {mustBeNumericOrLogical}  = false;        %--- Show the result
    options.SaveTheResult {mustBeNumericOrLogical}  = false;        %--- Save the result
    options.Plotsurfaces  {mustBeNumericOrLogical}  = false;        %--- Plot surfaces
end

%--- Material properties
constants.Tc      = x0(1);       % Curie temperature [K]
constants.gj      = x0(2);       % Landé factor
constants.J       = x0(3);       % Total angular momentum
constants.thetaD  = x0(4);       % Debye temperature (K)
constants.Ns      = x0(5);       % Number of spins pr unit mass in kg^-1
constants.M       = x0(6);       % Molar mass [kg/mol]
constants.gamma_e = x0(7);       % Sommerfeld constant (J/kgK^2)

%--- Magnetic flux density of Earth magnetic field
Bearth=0.00005;
if (min(Barr) < Bearth)
    disp('Minimum magnetic field is 0.00005 T')   
end

%--- Choose parameters to calculate [1 = calculate, 0 = do not calculate]
calculate_cp  = 1;
calculate_S   = 1;
calculate_mag = 1;
calculate_DT  = 1;

%------------------------ END OF USER INPUT ------------------------------------

%--- Physical constants
constants.Na      = 6.0221415e23;    % Avogadro's Number in 1/mol
constants.mu_0    = 4*pi*1e-7;       % Permeability of free space in N/A^2
constants.mu_b    = 9.2741e-24;      % Bohr magneton in J/T
constants.gamma_0 = 6.93e-2;         % Sommerfeld constant in J kg^-1 K^-2
constants.kb      = 1.3806505*1E-23; % Boltzmann's constant in J/K
    
%--- Derived constants
constants.N_atoms = constants.Na/constants.M;   % Number of atoms pr unit mass in kg^-1
   
%--- Define arrays
Cp      = zeros(length(Barr),length(Tarr));
DT      = Cp;
S       = Cp;
mag     = Cp;

%--- Lattice and electronic contribution to the total heat capacity
%--- These do not depend on constants.Tc and B and can thus be calculated outside the loop 
if (calculate_cp == 1)
    Cl = Clat(constants,Tarr);
    Ce = Cele(constants,Tarr);
end

%--- The mean field constant (which depends on constants.Tc)
constants.Nint=(3.*constants.kb.*constants.Tc)./(constants.Ns.*constants.gj.^2.*constants.mu_b.^2.*constants.J.*(constants.J+1)); 

%--- The zero field entropy
S0=Smag(constants,Bearth,Tarr)+Slat(constants,Tarr)+Sele(constants,Tarr);

for i = 1:length(Barr)
    if (calculate_cp == 1)
        %--- The specific heat capacity
        Cp(i,:) = Cmag(constants,Barr(i),Tarr)+Cl+Ce;
    end

    if (calculate_DT == 1)
        %--- The adiabatic temperature change
        DT(i,:) = DeltaTad(constants,S0,Barr(i),Tarr);
    end
  
    if (calculate_S == 1)
        %--- The specific entropy
        S(i,:)=Smag(constants,Barr(i),Tarr)+Slat(constants,Tarr)+Sele(constants,Tarr);
    end

    if (calculate_mag == 1)
        %--- The specific magnetization
        mag(i,:) = sigma(constants,Barr(i),Tarr);
    end
end

%--- The specific entropy change
DS = S-S0;

if (options.ShowTheResult)
    %--- Visualize the results
    if (calculate_cp == 1)
        [figureX,figX] = plot_results(Tarr,Cp);
        xlabel(figX,'Temperature [K]')
        ylabel(figX,'c_p [J kg^{-1} K^{-1}]')
        print('-dpng','Cp.png')
    end
    
    if (calculate_DT == 1)
        [figureX,figX] = plot_results(Tarr,DT(:,:));
        xlabel(figX,'Temperature [K]')
        ylabel(figX,'\Delta{}T_{ad} [K]')
        print('-dpng','DeltaTad.png')
    end
    
    if (calculate_S == 1)
        [figureX,figX] = plot_results(Tarr,DS(:,:));
        xlabel(figX,'Temperature [K]')
        ylabel(figX,'\Delta{}s_{iso} [J kg^{-1} K^{-1}]')
        print('-dpng','DeltaS.png')
        
        [figureX,figX] = plot_results(Tarr,S);
        xlabel(figX,'Temperature [K]')
        ylabel(figX,'S')
    end
    
    if (calculate_mag == 1)
        [figureX,figX] = plot_results(Tarr,mag);
        xlabel(figX,'Temperature [K]')
        ylabel(figX,'\sigma [emu/g or Am^2 kg^{-1}]')
        print('-dpng','Magnetization.png')
    end
end

if (options.Plotsurfaces)
    figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig1 = axes('Parent',figure1,'Layer','top','FontSize',16); axis off; grid off;
    surf(Tarr,Barr,Cp,'Linestyle','none')
    xlabel('Temperature [K]')
    ylabel('Magnetic flux density [T]')
    zlabel('Cp [J/(kg\cdot K)]')
    axis tight
    
    figure2= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig2 = axes('Parent',figure2,'Layer','top','FontSize',16); axis off; grid off;
    surf(Tarr,Barr,DT,'Linestyle','none')
    xlabel('Temperature [K]')
    ylabel('Magnetic flux density [T]')
    zlabel('\Delta T [K]')
    axis tight
    
    figure3= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig3 = axes('Parent',figure3,'Layer','top','FontSize',16); axis off; grid off;
    surf(Tarr,Barr,S,'Linestyle','none')
    xlabel('Temperature [K]')
    ylabel('Magnetic flux density [T]')
    zlabel('\Delta S')
    axis tight
    
    figure4= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto'); fig4 = axes('Parent',figure4,'Layer','top','FontSize',16); axis off; grid off;
    surf(Tarr,Barr,mag,'Linestyle','none')
    xlabel('Temperature [K]')
    ylabel('Magnetic flux density [T]')
    zlabel('Magnetization')
    axis tight
end

if (options.SaveTheResult)
    %-------------- Output the values as tables
    write_table_file(Tarr,Barr,Cp,'cp',strcat(['constants.Tc = ' num2str(constants.Tc) ]))
    write_table_file(Tarr,Barr,DT,'DeltaTad',strcat(['constants.Tc = ' num2str(constants.Tc) '\t Starting field = ' num2str(Bearth) ' T']))
    write_table_file(Tarr,Barr,S,'Entropy',strcat(['constants.Tc = ' num2str(constants.Tc) ]))
    write_table_file(Tarr,Barr,mag,'Magnetization',strcat(['constants.Tc = ' num2str(constants.Tc) ]))

    %--- Save the data in a matlab datafile
    save('MFT.mat', 'Tarr', 'Barr', 'Cp', 'DT', 'S', 'DS', 'mag', 'constants');
end

end

%% -----------------------------------------------------------------------
%% ------------------- The physical effects ------------------------------
%% -----------------------------------------------------------------------

%% Calculate the specific heat

%--- The specific electronic heat capasity [J/(kgK)]
function Ce = Cele(constants,T)
    Ce=constants.gamma_e.*T;
end

%--- The specific lattice heat capacity [J/(kgK)]
function Cl = Clat(constants,T)
    % Function to be integrated
    F = @(x) (x.^4.*exp(x))./((exp(x)-1).^2);
    
    % Loop over the temperature array
    for i=1:length(T)
        Cl(i)=9*constants.N_atoms*constants.kb*(T(i)./constants.thetaD).^3*quadl(F,1E-10,constants.thetaD./T(i));
    end
end

%--- The magnetic contribution to the heat capacity using the mean field theory [J/(kgK)]
function [Cm] = Cmag(constants,B,T)
    Cm=-B.*dsigmadT(constants,B,T)-0.5.*constants.Nint.*dsigma2dT(constants,B,T);
end

%% Calculate DeltaTad

%--- The adiabatic temperature change
%--- Calculating by using that the entropy in two fiels must be equal for an adiabatic process
function DT = DeltaTad(constants,S0,B2,T)
    DT = zeros(size(T));
    
    % Loop over the temperature array
    for i=1:length(T)
        F = @(x) (Smag(constants,B2,T(i)+x)+Slat(constants,T(i)+x)+Sele(constants,T(i)+x)-S0(i));
        DT(i)=fzero(F,[-1 100]);
    end
end

%% Calculate the specific entropy

%--- The magnetic entropy [J/(kgK)]
function Sm = Smag(constants,B,T)
     % Determination of the argument and solution to the Brillouin function
    [Bjchi,chi] = brillouinsolution(constants,B,T);

    % Calculation of the magnetic entropy
    logarg1=sinh(((2.*constants.J+1)./(2.*constants.J)).*chi);
    logarg2=sinh(chi./(2.*constants.J));
    Sm=constants.kb.*constants.Ns.*((log(logarg1)-log(logarg2))-chi.*Bjchi);
end

%--- The electronic entropy [J/(kgK)]
function Se = Sele(constants,T)
    Se=constants.gamma_e.*T;
end

%--- The lattice entropy [J/(kgK)]
function Sl = Slat(constants,T)
    % Function to be integrated
    G = @(x) x.^3./(exp(x)-1);

    Sl = zeros(size(T));
    
    % Loop over the temperature array
    for i=1:length(T)
        % Debye interpolation formula
        Sl(i)=constants.kb.*constants.N_atoms.*(-3.*log(1-exp(-constants.thetaD./T(i)))+12.*(T(i)./constants.thetaD).^3.*quadl(G,1E-10,constants.thetaD/T(i)));
    end
end

%% Calculate the specific magnetization
%--- The specific magnetization in the mean field theory [(A*m^2)/kg]
function sigma = sigma(constants,B,T)
    % Solution of the Brillouin function
    [Bjchi,chi] = brillouinsolution(constants,B,T);

    % Calculation of the specific magnetism
    sigma = constants.Ns*constants.gj*constants.J*constants.mu_b*Bjchi;
end

%% The Brillouin function and solution

% Function to determine the argument to and the result from the Brilloun
function [Bjchi,chi] = brillouinsolution(constants,B,T) % chi     : Argument to the Brilloun function, Bjchi   : Result from to the Brilloun function

    % Loop over the temperature array
    for i=1:length(T)
        % Residual function
        F = @(chi) (constants.gj.*constants.J.*constants.mu_b.*B)./(constants.kb.*T(i))+(3.*constants.Tc.*constants.J.*brillouin(constants,chi))./(T(i).*(constants.J+1))-chi;
        % Initializing bisection interval
        chi_low=-1;
        chi_high=1;

        % Determining high and low guess for numerical calculation of function roots
        while((sign(F(chi_low))==sign(F(chi_high))) & chi_high<1E6)
            chi_high=chi_high.*2;
        end

        chi(i)=fzero(F,[chi_low chi_high]);
        Bjchi(i)=brillouin(constants,chi(i));
    end
end

%--- The brillouin function
function Bjchi = brillouin(constants,chi)
    if chi<1E-8
        Bjchi=0;
    else
        Bjchi=((2.*constants.J+1)./(2.*constants.J)*coth(((2.*constants.J+1).*chi)./(2.*constants.J))-1./(2.*constants.J).*coth(chi./(2.*constants.J)));
    end
end

%% Calculate the derivative of the magnetization

%--- The first derivative of the specific magnetization with respect to the temperature [(A*m^2)/(kgK)]
function dsigmadT = dsigmadT(constants,B,T)
    % The first derivative is determined numerically 
    h=1E-1;
    dsigmadT=(sigma(constants,B,T+h)-sigma(constants,B,T-h))./(2*h);
end

%--- The second derivative of the specific magnetization squared with respect to temperature  [(A*m^2)/(kgK)]
function d2sigmadT = dsigma2dT(constants,B,T)
    % The second derivative is determined numerically 
    h=1E-1;
    d2sigmadT = (sigma(constants,B,T+h).^2-sigma(constants,B,T-h).^2)./(2*h);
end


%% -----------------------------------------------------------------------
%% ------------------- Auxillary functions -------------------------------
%% -----------------------------------------------------------------------

%% Visualize the data
function [figure1,fig1] = plot_results(Xarr,Yarr)
    figure1 = figure('PaperSize',[20.98 29.68]);
    fig1 = axes('Parent',figure1,'Layer','top','FontSize',18);
    grid on
    hold all
    
    Marker_arr = {'x' 's' 'd' 'o' '<' '>' 'v' '^'};
    
    for i = 1:length(Yarr(:,1))
        plot(Xarr,Yarr(i,:),'Marker',Marker_arr{i},'Markersize',8,'Linestyle','none','color','k');
    end
    
    xlim([min(Xarr) max(Xarr)]);
end

%% Write table with data to file
function write_table_file(Temp,H_field,data_z,data_type,header)
    %--- Write the result to file
    fid = fopen(strcat(['MFT_' data_type '_table.txt']),'w');
    
    for p=1:length(data_z(:,1))
        %--- First the header and a row containing the temperature
        if (p == 1)
            fprintf(fid,strcat(['%% ' header '\r\n']));
            fprintf(fid,strcat(['%%H field first colomn\tTemp row\t' data_type ' data array\r\n']));
            for j = 1:length(Temp)
                if (j == 1)
                    fprintf(fid,'%6.4f\t\t',0.000);   %--- Print a zero to begin with
                end
                fprintf(fid,'%6.4f\t\t',Temp(j));
                if j == length(Temp)
                    fprintf(fid,'\r\n');
                end
            end
        end
        
        %--- Print the field
        fprintf(fid,'%6.4f\t\t',H_field(p));
        
        %--- Print the data
        for j = 1:length(data_z(p,:))
            fprintf(fid,'%8.4f\t\t',data_z(p,j));
            if j == length(data_z(p,:))
                fprintf(fid,'\r\n');
            end
        end
    end
    fclose(fid);
end
