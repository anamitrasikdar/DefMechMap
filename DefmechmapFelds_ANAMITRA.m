tic;
% close all;
% clear all;
%%  Constants' input section

R = 8.31446;        % Gas constant
F = 534;           % H2O fugacity 
r = 1;                  % Fugacity exponent

% Temperature inpute subsection
prompt = {'Enter the temperature in °C:'}; 
dlg_title = 'Temperature input';
num_lines = 1;
default = {'400'};             % default value for temperature
t = str2double(inputdlg(prompt, dlg_title, num_lines, default)); % opening dialog box
T = t+273;                      % conversion from degree centrigrade to Kelvin
P = 700;                          %Mpa
% %Feldspar creep parameters (2004 Rybacki and Dresen)
% % disl => dislocation related constants
% ndisl = 3;                       % stress exponent
% mdisl = 0;                      % grain size exponent
% Qdisl = 332000;             %J mol^-1 % internal energy constant
% Adisl = 2.511e3;             %MPa^-n mu^m s^-1 % material constant
% Cdisl = exp((-Qdisl)/(R*T)); % derived constant
% ACFdisl = Adisl*Cdisl*(F^r); % derived constant

%Feldspar creep parameters (2006 Rybacki et al JGRSE)
% disl => dislocation related constants
V = 38;                    % cm3/mol
ndisl = 3;                       % stress exponent
mdisl = 0;                      % grain size exponent
Qdisl = 345000;             %J mol^-1 % internal energy constant
Adisl = 10^0.2;             %MPa^-n mu^m s^-1 % material constant
Cdisl = exp((-Qdisl-P*V)/(R*T)); % derived constant
ACFdisl = Adisl*Cdisl*(F^r); % derived constant


%Feldspar diffusion parameters (2006 Rybacki et al JGRSE)
% dif => diffusion related constants
ndif = 1;                        % stress exponent
mdif = 3;                       % grain size exponent
Qdif = 159000;               % J mol^-1 % internal energy constant
Adif = 10^-0.7;             % MPa^-n mu^m s^-1 % material constant
Cdif = exp((-Qdif-P*V)/(R*T)); % derived constant
ACFdif = Adif*Cdif*(F^r); % derived constant

% %Feldspar diffusion parameters (2004 Rybacki and Dresen)
% % dif => diffusion related constants
% ndif = 1;                        % stress exponent
% mdif = 3;                       % grain size exponent
% Qdif = 193000;               % J mol^-1 % internal energy constant
% Adif = 7.9433e3;             % MPa^-n mu^m s^-1 % material constant
% Cdif = exp((-Qdif)/(R*T)); % derived constant
% ACFdif = Adif*Cdif*(F^r); % derived constant

% Cardano correction factor for differential stress and strain rate
stressCardano = nthroot(3,2);     % sigma = nthroot(3,2)*tau; 
strainCardano = 1/(nthroot(3,2)); % gamma_dot = nthroot(3,2)*eplison_dot;

% Boundary constants for defining the diffusion and dislocation regime 
BC = ACFdisl/ACFdif;

%% Calculation section 

% Calculations of deformation mechanism map

d = 100 : 100000;          % range effectively [1:1000 µm], grain size reduced for fine discretization
stress = 1 : 1000;          % range effectively [1:1000 MPa], stress not discretised for insufficient memory
felstrainrdisl = [length(stress), length(d)];
felstrainrdif =  [length(stress), length(d)];
felstrainrtot=  [length(stress), length(d)];
%  bar = waitbar(100, 'Program is runnung, please wait' );
for j = 100 : length(d)

    for i = 1 : length(stress)
        felstrainrdisl(i, j) = ACFdisl*((i*stressCardano)^ndisl)*((j/100)^-mdisl);
        felstrainrdif(i, j) = ACFdif*((i*stressCardano)^ndif)*((j/100)^-mdif);
        felstrainrtot(i, j) =(felstrainrdisl(i, j)+felstrainrdif(i, j))/strainCardano;
     
    end
     
end

% Calculations of boundary(s) of deformation mechanism map
% CONDITION => epsilon_disl = epsilon_dif
    for stress = 1 : 1000           %range effectively [1:1000]
        rexd = [1, length(stress)];
        diffstress = [1, length(stress)];
        rexd (stress) = ((BC*(stress*stressCardano)^(ndisl-ndif))^(-1/(mdif-mdisl)))*100;
        diffstress (stress) = stress*stressCardano;
    
    end 
    fprintf('DMM and Boundary - Calculated\n'); % optional
%% Optional error checking for boundary(s)
    felstrainrdif1 = ACFdif*((diffstress).^ndif).*((rexd/100).^-mdif);
    felstrainrdisl1 = ACFdisl*((diffstress).^ndisl).*((rexd/100).^-mdisl);
    test = felstrainrdif1./felstrainrdisl1;

%% Plot section
% Plotting of feldspar deformation mechanism map
    fprintf('Plotting started XO XO\n'); % optional
figure
set (gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % fullpage view
h=contour(felstrainrtot,[1e-10 1e-11 1e-12 1e-13 1e-14 1e-15],'-k','ShowText','off'); % countour of strain rate
clabel (h);                                             % contour handle control
hold on
plot (rexd,diffstress./stressCardano,'-g'); % boundary of dislocation and diffusion
set(gca,'xscale','log')
set(gca,'yscale','log')                              % set axes for log scale
title(sprintf('Feldspar DMM at T =  %i °C', T-273)) %conversion back to °C
xlabel('d (\mum) / 100');
ylabel('\sigma (MPa)');
 xlim ([1e2 1e5])                                    %xaxis in reality [x value/100 µm]
 ylim ([1e1 1e3])                                    %yaxis in reality [MPa]

%%
% Data input for scatter plot of sample points
diff_stress = [18.7, 16.5, 17.6];
d_rms = [20.2, 23, 21.5]*100;
col_lab = [1 0 0; 0 1 0; 0 0 1];
scatter(d_rms, diff_stress, 100, col_lab, 'filled', 'DisplayName', 'SB2')

diff_stress2 = [17.3, 23.4, 16.6];
d_rms2 = [21.9, 16.1, 22.7].*100;
col_lab2 = [1 0 0;  0 1 0;  0 0 1];
scatter(d_rms2, diff_stress2, 100, col_lab2, 's', 'DisplayName', 'SB1')
hold off
toc;

%% Export file to desired format

fig_path = 'F:\Anamitra\Work\Deformation Mechanism Map\';
sample_name = ('ALL');
phase_name = ('Feldspar');
filetype1 = ('.pdf');
filetype2 = ('.png');
fig_name1 = strcat(phase_name,'_',sample_name,'_', num2str(t), filetype1);
fig_name2 = strcat(phase_name,'_',sample_name,'_', num2str(t), filetype2);
filename1 = fullfile(fig_path, fig_name1);
filename2 = fullfile(fig_path, fig_name2);
exportgraphics(gcf,filename2,'Resolution',300);
exportgraphics (gcf , filename1, 'ContentType', 'vector');