tic;
close all;
clearvars;

%%  Constants' input section

% Input subsection
prompt = {'Enter pressure in MPa:','Enter the temperature in °C:', 'Enter the fugacity factor'}; 
dlg_title = 'Input P-T and fug_factor';
dims = [1 70];
default = {'450', '400', '0.229'};  % Default values
answer = inputdlg(prompt, dlg_title, dims, default);
P = str2double(answer{1});  %Mpa
T = str2double(answer{2}) + 273;    % Conversion from degree centrigrade to Kelvin
f = str2double(answer{3}); % Fugacity factor

% Other constants in the flow laws    
R = 8.31446;    % Gas constant
F = P*f;    % H2O fugacity 
r = 1;  % Fugacity exponent

% Quartz creep parameters (2001 HIrth et. al)
% disl => dislocation related constants
V = 0;                            % cm3/mol    
ndisl = 4;                       % Stress exponent
mdisl = 0;                      % Grain size exponent
Qdisl = 135000;              % J mol^-1 % internal energy constant
Adisl =  6.3*10^-(12);      % MPa^-n mu^m s^-1 % material constant
Cdisl = exp((-Qdisl-P*V)/(R*T)); % Derived constant
ACFdisl = Adisl*Cdisl*(F^r); % Derived constant

%Quartz diffusion parameters (2000 Brodie and Rutter)
% dif => diffusion related constants
ndif = 1;                        % Stress exponent
mdif = 2;                       % Grain size exponent
Qdif = 220000;               % kJ mol^-1 % internal energy constant
Adif = 10^(-0.2);             % MPa^-n mu^m s^-1 % material constant
Cdif = exp((-Qdif-P*V)/(R*T)); % Derived constant
ACFdif = Adif*Cdif*(F^r); % Derived constant

% Cardano correction factor for differential stress and strain rate
stressCardano = nthroot(3,2);     % sigma = nthroot(3,2)*tau; 
strainCardano = 1/(nthroot(3,2)); % gamma_dot = nthroot(3,2)*eplison_dot;

% Boundary constants for defining the diffusion and dislocation regime 
BC = ACFdisl/ACFdif;

%% Calculation section 

% Calculations of deformation mechanism map

d = 10 : 100000;          % range effectively [1:1000 µm], grain size reduced for fine discretization
stress = 1 : 1000;          % range effectively [1:1000 MPa], stress not discretised for insufficient memory
qtzstrainrdisl = [length(stress), length(d)];   % Intiating variables
qtzstrainrdif =  [length(stress), length(d)];   % Intiating variables
qtzstrainrtot=  [length(stress), length(d)];    % Intiating variables

nLoops = length(d); % WAITBAR (OPTIONAL)
progressbarText(0, 50); % WAITBAR (OPTIONAL)

for j = d(1) : length(d)

loopCnt =j;  % WAITBAR (OPTIONAL)
 
    for i = stress(1) : length(stress)
        qtzstrainrdisl(i, j) = ACFdisl*((i*stressCardano)^ndisl)*((j/100)^-mdisl);
        qtzstrainrdif(i, j) = ACFdif*((i*stressCardano)^ndif)*((j/100)^-mdif);
        qtzstrainrtot(i, j) =(qtzstrainrdisl(i, j)+qtzstrainrdif(i, j))/strainCardano;
     
    end
    
progressbarText(loopCnt/nLoops);    % WAITBAR (OPTIONAL)
  
end

% Calculations of boundary(s) of deformation mechanism map
% CONDITION => epsilon_disl = epsilon_dif
 rexd = zeros(1, length(stress));   % Initiating variables
        diffstress = rexd;                % Initiating variables  
    for stress = 1 : 1000           %range effectively [1:1000]
        rexd (stress) = ((BC*(stress*stressCardano)^(ndisl-ndif))^(-1/(mdif-mdisl)))*100;
        diffstress (stress) = stress*stressCardano;
    
    end 
    fprintf('DMM and Boundary - Calculated\n'); % Optional
    
%% Optional error checking for boundary(s)
%     qtzstrainrdif1 = ACFdif*((diffstress).^ndif).*((rexd/100).^-mdif);
%     qtzstrainrdisl1 = ACFdisl*((diffstress).^ndisl).*((rexd/100).^-mdisl);
%     test = qtzstrainrdif1./qtzstrainrdisl1;

%========================
% test = [1 1 1 1 1 ... 1]    check the test variable
%========================

%% Plot section
% Plotting of feldspar deformation mechanism map
  
fprintf('Plotting started...  (O_o) ( .•_•)\n'); % Optional

figure
% Countour of strain rate

range = [1e-13 1e-12 1e-11 1e-10 1e-9 1e-8]; % Select the strain contours
h=contour(qtzstrainrtot, range, 'ShowText','off', 'LineWidth', 2);  % Contour plot
clabel(h);                                             % Contour label
colormap(copper);  
set(gca, 'ColorScale', 'log'); % Assign contour colormap, Otherwise replace the following line
% h=contour(qtzstrainrtot, range, '-k' 'ShowText','off', 'LineWidth', 2);
% hcb = colorbar;
% hcb.Ruler.Exponent =10;

set(gca, 'clim', [min(range) max(range)])
set (gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);  % fullpage view

hold on
plot (rexd,diffstress./stressCardano,'-m'); % boundary of dislocation and diffusion

set(gca,'xscale','log')
set(gca,'yscale','log')                              % set axes for log scale


title(sprintf('Qtz DMM at P = %iMPa and T =  %i°C', P, T-273)) %conversion back to °C
xlabel('d (\mum) / 100'); % Apply effective range
ylabel('\sigma (MPa)');
 
xlim ([1e2 1e5])                                    %xaxis in reality [x value/100 µm]
ylim ([1e1 1e3])                                    %yaxis in reality [MPa]

%% Input the diff stress vs grain size data points from the EBSD analysis
% Data input for scatter plot of sample points
col_lab = [1 0 0; 0 1 0; 0 0 1]; % R, G, B=> T, B, M respectively

diff_stress = [70, 68.4, 73.1];   % MPa
d_rms = [20.4, 21, 19.1].*100;  % Adjust for the effective range of grain size
scatter(d_rms, diff_stress, 100, col_lab, 's',  'DisplayName', 'SB1')

diff_stress2 = [76.5, 81, 63.7];      % MPa
d_rms2 = [17.9, 16.5, 23.2].*100;  % Adjust for the effective range of grain size
scatter(d_rms2, diff_stress2, 100, col_lab, 'filled', 'DisplayName', 'SB2')

hold off
toc;
