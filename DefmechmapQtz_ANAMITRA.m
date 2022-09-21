tic;
% close all;
% clear all;
%%  Constants' input section

R = 8.31446;        % Gas constant
F = 209;           % H2O fugacity 
r = 1;                  % Fugacity exponent

% Temperature inpute subsection
prompt = {'Enter the temperature in °C:'}; 
dlg_title = 'Temperature input';
num_lines = 1;
default = {'400'};             % default value for temperature
t = str2double(inputdlg(prompt, dlg_title, num_lines, default)); % opening dialog box
T = t+273;                      % conversion from degree centrigrade to Kelvin
P = 0;
% Quartz creep parameters (2001 HIrth et. al)
% disl => dislocation related constants
V = 0;                       %cm3/mol    
ndisl = 4;                       % stress exponent
mdisl = 0;                      % grain size exponent
Qdisl = 135000;             %J mol^-1 % internal energy constant
Adisl =  6.3*10^-(12);           %MPa^-n mu^m s^-1 % material constant
Cdisl = exp((-Qdisl-P*V)/(R*T)); % derived constant
ACFdisl = Adisl*Cdisl*(F^r); % derived constant

% % Quartz creep parameters (2019 Tokle et al)
% % disl => dislocation related constants
% V = 0;                       %cm3/mol    
% ndisl = 4;                       % stress exponent
% mdisl = 0;                      % grain size exponent
% Qdisl = 140000;             %J mol^-1 % internal energy constant
% Adisl = 8*10^(-12);           %MPa^-n mu^m s^-1 % material constant
% Cdisl = exp((-Qdisl-P*V)/(R*T)); % derived constant
% ACFdisl = Adisl*Cdisl*(F^r); % derived constant

% % Quartz creep parameters (2018 Fukuda et al)
% % disl => dislocation related constants
% V = 0;                       %cm3/mol    
% ndisl = 1.7;                       % stress exponent
% mdisl = 0.51;                      % grain size exponent
% Qdisl = 183000;             %J mol^-1 % internal energy constant
% Adisl = 10^(-2.97);           %MPa^-n mu^m s^-1 % material constant
% Cdisl = exp((-Qdisl-P*V)/(R*T)); % derived constant
% ACFdisl = Adisl*Cdisl*(F^r); % derived constant

%Quartz diffusion parameters (2000 Brodie and Rutter)
% dif => diffusion related constants
ndif = 1;                        % stress exponent
mdif = 2;                       % grain size exponent
Qdif = 220000;               % kJ mol^-1 % internal energy constant
Adif = 10^(-0.2);             % MPa^-n mu^m s^-1 % material constant
Cdif = exp((-Qdif-P*V)/(R*T)); % derived constant
ACFdif = Adif*Cdif*(F^r); % derived constant

% Cardano correction factor for differential stress and strain rate
stressCardano = nthroot(3,2);     % sigma = nthroot(3,2)*tau; 
strainCardano = 1/(nthroot(3,2)); % gamma_dot = nthroot(3,2)*eplison_dot;

% Boundary constants for defining the diffusion and dislocation regime 
BC = ACFdisl/ACFdif;

%% Calculation section 

% Calculations of deformation mechanism map

d = 10 : 100000;          % range effectively [1:1000 µm], grain size reduced for fine discretization
stress = 1 : 1000;          % range effectively [1:1000 MPa], stress not discretised for insufficient memory
qtzstrainrdisl = [length(stress), length(d)];
qtzstrainrdif =  [length(stress), length(d)];
qtzstrainrtot=  [length(stress), length(d)];

for j = d(1) : length(d)

    for i = stress(1) : length(stress)
        qtzstrainrdisl(i, j) = ACFdisl*((i*stressCardano)^ndisl)*((j/100)^-mdisl);
        qtzstrainrdif(i, j) = ACFdif*((i*stressCardano)^ndif)*((j/100)^-mdif);
        qtzstrainrtot(i, j) =(qtzstrainrdisl(i, j)+qtzstrainrdif(i, j))/strainCardano;
     
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
%     qtzstrainrdif1 = ACFdif*((diffstress).^ndif).*((rexd/100).^-mdif);
%     qtzstrainrdisl1 = ACFdisl*((diffstress).^ndisl).*((rexd/100).^-mdisl);
%     test = qtzstrainrdif1./qtzstrainrdisl1;

%% Plot section
% Plotting of feldspar deformation mechanism map
    fprintf('Plotting started XO XO\n'); % optional
figure
set (gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % fullpage view
h=contour(qtzstrainrtot,[1e-8 1e-9 1e-10 1e-11 1e-12 1e-13],'-k','ShowText','off'); % countour of strain rate
clabel (h);                                             % contour handle control
hold on
plot (rexd,diffstress./stressCardano,'-g'); % boundary of dislocation and diffusion
set(gca,'xscale','log')
set(gca,'yscale','log')                              % set axes for log scale
title(sprintf('Quartz DMM at T =  %i °C', T-273)) %conversion back to °C
xlabel('d (\mum) / 100');
ylabel('\sigma (MPa)');
 xlim ([1e1 1e5])                                    %xaxis in reality [x value/100 µm]
 ylim ([1e1 1e3])                                    %yaxis in reality [MPa]

%%
%Data input for scatter plot of sample points
diff_stress = [76.5, 81, 63.7];
d_rms = [17.9, 16.5, 23.2].*100;
col_lab = [1 0 0; 0 1 0; 0 0 1];
scatter(d_rms, diff_stress, 100, col_lab, 'filled', 'DisplayName', 'SB2')

diff_stress2 = [70, 68.4, 73.1];
d_rms2 = [20.4, 21, 19.1].*100;
col_lab2 = [1 0 0;  0 1 0;  0 0 1];
scatter(d_rms2, diff_stress2, 100, col_lab2, 's',  'DisplayName', 'SB1')
hold off
toc;

%% Export file to desired format

fig_path = 'F:\Anamitra\Work\Deformation Mechanism Map\';
sample_name = ('ALL');
phase_name = ('Quartz');
filetype1 = ('.pdf');
filetype2 = ('.png');
fig_name1 = strcat(phase_name,'_',sample_name,'_', num2str(t), filetype1);
fig_name2 = strcat(phase_name,'_',sample_name,'_', num2str(t), filetype2);
filename1 = fullfile(fig_path, fig_name1);
filename2 = fullfile(fig_path, fig_name2);
exportgraphics(gcf,filename2,'Resolution',300);
exportgraphics (gcf , filename1, 'ContentType', 'vector');