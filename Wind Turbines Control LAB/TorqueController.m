clear all
close all
clc

%% Load Performance Data
load('Performance_data_15MW.mat')

%% Initialize Parameters
R = 120;                                    % Rotor radius (meters)
rho_air = 1.225;                            % Air density (kg/m^3)
tau_gb = 1;                                 % Gearbox ratio
eta_gb = 1;                                 % Generator efficiency
U = 10:0.1:11;                              % Wind speed range (m/s)
TSR = Rotor_Lamda';                         % Tip-speed ratio
betha = linspace(-5, 5, size(Rotor_cP, 1)); % Blade Pitch angle

% Meshgrid for plotting
[x, y] = meshgrid(TSR, betha);
surf(x, y, Rotor_cP)
xlabel('Tip-speed ratio TSR [-]')
ylabel('Blade pitch')
zlabel('Power coefficient [-]')
grid on

%% Find the Maximum Power Coefficient
[Cp_max, linearIdx] = max(Rotor_cP(:));

[index_betha, index_TSR] = ind2sub(size(Rotor_cP), linearIdx);
betha_opt = y(index_betha, 1);
TSR_opt = x(1, index_TSR);
fprintf('The max power coefficient Cp is %s\n', num2str(Cp_max))

%% Plot the Power Coefficient Contour
figure 
contourf(x, y, Rotor_cP)
colorbar('eastoutside');
c = colorbar;  
c.Label.String = 'Power coefficient [-]'; 
c.Label.FontSize = 11;
c.Label.FontWeight = 'bold';
title('Fine pitch and optimal TSR');
xlabel('Tip-Speed Ratio TSR [-]')
ylabel('Blade pitch')
hold on
plot(x(1, index_TSR), y(index_betha, 1), '*r', 'MarkerSize', 10, 'DisplayName', 'Cp max');
legend('Location', 'best');

%% Rated Power and Rotor Speed
figure
rated_power = 15e6;  % 15 MW rated power
for ii = 1:length(U)
    Rotor_speed = TSR .* U(ii) / R * (60 / (2 * pi));   % Rotor speed in rpm
    Pa = 0.5 * rho_air * U(ii)^3 * pi * R^2 * Rotor_cP(index_betha, :);  % Aerodynamic power
    
    % Plot aerodynamic power
    plot(Rotor_speed, Pa * 1e-6, 'DisplayName', sprintf('U = %.2f m/s', U(ii))); 
    hold on; grid on;
    
    % Check if the power exceeds rated 15 MW
    [max_P, idx_max_P] = max(Pa);
    if max_P >= rated_power && ~exist('rated_wind_speed', 'var')
        rated_wind_speed = U(ii);
        rated_omega = Rotor_speed(idx_max_P);  % Rotor speed at rated power
        rated_P = max_P;  % Rated aerodynamic power
    end
end

% Add dashed lines for rated values
xline(rated_omega, '--r', 'LineWidth', 1, 'DisplayName', 'Rated Rotor Speed');
yline(rated_P * 1e-6, '--r', 'LineWidth', 1, 'DisplayName', 'Rated Power');

% Labels and legend
xlabel('Rotor speed [rpm]')
ylabel('Aerodynamic rotor power [MW]')
legend('Location', 'northeastoutside');
title('Rated Power and Rotor Speed')
fprintf('Rated power is: %.2f MW at rotor speed: %.2f rpm and wind speed: %.2f m/s\n', rated_P * 1e-6, rated_omega, rated_wind_speed);

%% Rated torque and Rotor Speed
figure
for ii = 1:length(U)
    Rotor_speed = TSR .* U(ii) / R * (60 / (2 * pi));                    % Rotor speed in rpm
    Pa = 0.5 * rho_air * U(ii)^3 * pi * R^2 * Rotor_cP(index_betha, :);  % Aerodynamic power
    Generator_speed = Rotor_speed.*tau_gb;                               % Generator speed for the x axis
    Pg = Pa*eta_gb;                                                      % Power at the generator   
    Qg = Pa./(Generator_speed*2*pi/60);                                  % Generator Torque wrt generator speed in rad/s

    % Plot generator torque
    plot(Generator_speed, Qg * 1e-3,'k','DisplayName', sprintf('U = %.2f m/s', U(ii))); 
    hold on; grid on;
    xlim([0 8])
end


% rated generator speed and rated generator torque
rated_omega_gb_ideal = rated_omega*tau_gb;
rated_omega_gb_real = rated_omega_gb_ideal * 1.02; % real rated generator speed w gain of 2%
rated_Qb_ideal = rated_P * eta_gb / (rated_omega_gb_ideal*2*pi/60);
rated_Qb_real = rated_P * eta_gb / (rated_omega_gb_real*2*pi/60); 

% Add dashed lines for rated values
xline(rated_omega_gb_ideal, '--b', 'LineWidth', 1, 'DisplayName', 'Rated Generator Speed (ideal)');
yline(rated_Qb_ideal * 1e-3, '--b', 'LineWidth', 1, 'DisplayName', 'Rated Generator Torque (ideal)');
xline(rated_omega_gb_real, '--r', 'LineWidth', 1, 'DisplayName', 'Rated Generator Speed (real)');
yline(rated_Qb_real * 1e-3, '--r', 'LineWidth', 1, 'DisplayName', 'Rated Generator Torque (real)');

% Labels and legend
xlabel('Generator speed [rpm]')
ylabel('Generator Torque [KNm]')
legend('Location','northeastoutside');
title('Generator Torque and Genrator Speed - Optimal Gain')

% Plot the optimal generator gain
hold on
k_opt = 0.5 * rho_air * pi * R^5 * Cp_max / TSR_opt^3 / tau_gb^3 * eta_gb;  % Optimal generator gain
Qg_opt = k_opt * (Generator_speed*2*pi/60).^2; % Generator torque curve for maximum efficiency
plot(Generator_speed, Qg_opt * 1e-3, 'k','LineWidth',1.5,'DisplayName', 'Optimal Generator Curve');

%% Define transition speed from FAST tool
omega_A = 5;            % [Rpm] minimum generator speed
omega_B = 5.2;          % [Rpm] transition speed 1
omega_B2 = 6.5;         % [Rpm] transition speed 2
omega_B2 = 7.4;         % [Rpm] transition speed 2 --> enlarged optimal region
omega_C = 7.6;          % [Rpm] maximum generator speed

xline(omega_A, '--m', 'LineWidth', 1, 'DisplayName', 'Minimum Gen. speed');
xline(omega_B, '--m', 'LineWidth', 1, 'DisplayName', 'Transition speed 1');
xline(omega_B2, '--m', 'LineWidth', 1, 'DisplayName', 'Transition speed 2');
xline(omega_C, '--m', 'LineWidth', 1, 'DisplayName', 'Maximum Gen. speed');

%% Indexes of the speeds and coefficient
A=find(min(abs(omega_A-Generator_speed))==abs(omega_A-Generator_speed));
B=find(min(abs(omega_B-Generator_speed))==abs(omega_B-Generator_speed));
B2=find(min(abs(omega_B2-Generator_speed))==abs(omega_B2-Generator_speed));
C=find(min(abs(omega_C-Generator_speed))==abs(omega_C-Generator_speed));

S1=(Qg_opt(B)-Qg_opt(A))/(Generator_speed(B)-Generator_speed(A))/2;
S2=(rated_Qb_real-Qg_opt(B2))/(Generator_speed(C)-Generator_speed(B2+1));

%% Piecewise Function   
Qg_real(1:A)=0;
for ii=A+1:B
Qg_real(ii)=S1*(Generator_speed(B)-Generator_speed(A)/(B-A))*(ii-A);
end
Qg_real(B+1:B2)=Qg_opt(B+1:B2);

for ii=B2+1:C
    Qg_real(ii)=Qg_opt(B2)+S2*(Generator_speed(ii)-Generator_speed(B2));
end

 Qg_real(C+1:201)=0;
 for ii=C+1:201
   Pa = 0.5 * rho_air * U(5)^3 * pi * R^2 * Rotor_cP(index_betha, :);
   Qg_real(ii)=Pa(ii)/(Generator_speed(ii)*2*pi/60);
 end
hold on
plot(Generator_speed,Qg_real*1e-3,'LineWidth',2.5, 'DisplayName', 'Optimal Mode gain','Color',"c")
xregion(omega_A,omega_B,'FaceColor','y','DisplayName','Initial transition speed range')
xregion(omega_B,omega_B2,'FaceColor','g','DisplayName','Optimal speed range')
xregion(omega_B2,omega_C,'FaceColor','r','DisplayName','End transition speed range')

%% S1 and S2 for simulink
S1_sim = (k_opt*(omega_B*2*pi/60)^2 - 0)/(omega_B - omega_A)
S2_sim = (rated_Qb_real-k_opt*(omega_B2*2*pi/60)^2) / (omega_C - omega_B2)
%k_opt*(omega_B2*2*pi/60)^2 + S2_sim * (omega_C-omega_B2);