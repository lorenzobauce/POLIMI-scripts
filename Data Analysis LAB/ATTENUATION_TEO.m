clear; close all; clc

%% LOADING TF FROM SC EMS WITHOUT Rsh CIRCUIT
% visual check on the estimated transfer function
load("transfer_fun.mat")
H_sc=G_jk;
figure
title('Experimental Trasfer Function')
semilogy(abs(G_jk), 'LineWidth',1.5)
grid on
axis tight

% find f nat indexes
[~, locus_hsc] = findpeaks(log(abs(H_sc)),'MinPeakProminence',2);

% data coming from the experiments in short and open circuits
% nomenclature reffered to the paper
omega_itilde=[ 174.46  895.98  2567.5];
omega_i=[ 169.65  892.63  2554.3];

%% Piezoelectric path capacitance
% piezoelectric patch capacitance of tne blocked structure (i.e in
% sismographic region) evaluated by the analytical exact form
fig = openfig('cpi.fig');
fig = gcf;
ax = gca;
lines = findobj(ax, 'Type', 'line');
for i = 1:length(lines)
    frequency = get(lines(i), 'XData');
    capacity = get(lines(i), 'YData');
end

[~,index]=findpeaks(capacity,'MinPeakDistance', 100,'MinPeakHeight',4e-8);
mean_point1=floor((index(1)+index(2))/2);
mean_point2=floor((index(2)+index(3))/2);
mean_point3=floor((index(3)+index(4))/2);

Cp_1=capacity(mean_point1);
Cp_2=capacity(mean_point2);
Cp_3=capacity(mean_point3);

%% Building the shunt impedance
%Cp_1 = 37.5e-9; % [Farad] value used for our experiment

Beta2 = 0.7; 
C2_serie = Cp_1/Beta2; 

k1 = sqrt((omega_itilde(1)^2 - omega_i(1)^2)/omega_i(1)^2);
omega_sc = omega_i(1)*sqrt(1-(Beta2*k1^2 / (1-Beta2)));
omega_oc = omega_itilde(1);
omega_f = sqrt((omega_oc.^2 + omega_sc.^2)./2);
tau_opt = 1./omega_f;

Rsh = tau_opt * (1 - Beta2) / Cp_1;
EMEMCF1 = k1/sqrt((1-Beta2));
X1 = EMEMCF1 * omega_i * sqrt(Cp_1);

%% Compare the experimenatal results with the theoretical analytical effect coming from the single mode
% Performance evaluation of single mode control [mode 1]
% this function is the FRF evaluated with the optimal tau constant, function of:
% omega_i --> natural frequency of the mechanical structure
% omega_s
% omega_oc
% xii     --> structural damping of the structure 
% tau_e   --> electrical time constant which depends on the mode considered [mode 1]

Ceq = Cp_1*C2_serie/(C2_serie-Cp_1); % Equivalent capacitance for series config, evaluated for Cp in 1st mode

omega_step =(100:1/3:250);
f_step = 1/30;
o_step = 2*pi*f_step;
f_start = round(100/o_step/pi);
f_fine = round(250/o_step/pi);

xii = 0.0145; % real value evaluated from the first experiments in short circuit
H1 = @(Omega) (1+1i*tau_opt(1).*Omega)./(omega_sc(1).^2 - (1 + 2*xii*omega_i(1)*tau_opt(1)).*Omega.^2 + 1i.*Omega.*(tau_opt(1)*omega_oc(1).^2 + 2*xii*omega_i(1) - tau_opt(1).*Omega.^2));

figure
%title('FRF NC layout w optimal Rsh on 1^{st} mode H1(\Omega,\xi)')
%semilogy(100:o_step*pi:250, abs(H_sc(f_start:f_fine-1)),'k--', 'LineWidth',1)
hold on
grid on
semilogy(omega_step, abs(H1(omega_step)),'--','LineWidth',1.5);

for xii = [0.0001 0.001 0.01 0.1]
    H1 = @(Omega) (1+1i*tau_opt(1).*Omega)./(omega_sc(1).^2 - (1 + 2*xii*omega_i(1)*tau_opt(1)).*Omega.^2 + 1i.*Omega.*(tau_opt(1)*omega_oc(1).^2 + 2*xii*omega_i(1) - tau_opt(1).*Omega.^2));
    semilogy(omega_step, abs(H1(omega_step)), 'LineWidth',1); 
end

title('FRFs comparison on 1^{st} resonance frequency')
ylabel('Magnitude | |')
xlabel('Frequency [rad/s]')
legend('FRF NC series R_{sh} \xi estimated = 1.45%','FRF NC series R_{sh} \xi = 0.01%', 'FRF NC series R_{sh} \xi = 0.1%', 'FRF NC series R_{sh} \xi = 1%', 'FRF NC series R_{sh} \xi = 10%',Location='best')
hold off

%% Comparison teo
xii = 0.01;
Hsc = 1/(2*xii*omega_i(1)^2*sqrt(1-xii^2));
H1 = @(Omega) (1+1i*tau_opt(1).*Omega)./(omega_sc(1).^2 - (1 + 2*xii*omega_i(1)*tau_opt(1)).*Omega.^2 + 1i.*Omega.*(tau_opt(1)*omega_oc(1).^2 + 2*xii*omega_i(1) - tau_opt(1).*Omega.^2));
Adb_1 = 20*log10(Hsc./abs(H1(omega_f(1))));
disp('Attenuazione teorica calcolato come rapporto tra'); 
disp(['FRF elemento EMS in short circuit e con Rsh teorica: ', num2str(Adb_1), 'dB']);

%% Alternative formulation
K = k1*sqrt(Beta2/(1-Beta2));
Adb_teo_1 = 20*log10((EMEMCF1^2+2*sqrt(2)*xii*sqrt(2+EMEMCF1^2-2*K^2))/(4*xii*sqrt(1-xii^2)));
disp(['Attenuazione con formula teorica: ', num2str(Adb_teo_1), 'dB']);

%% Comparison exp
load("close.mat")
Hsc = abs(G_jk(locs(1)));
load("circuit.mat")
Hsh = abs(G_jk(locs(1)));
Adb_exp = 20*log10(Hsc/Hsh);
disp(['Attenuazione sperimentale: ', num2str(Adb_exp),'dB']);

