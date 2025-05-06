% Linearized aerodynamic model of the IEA 15MW 
% clear
close all
clc

%% Load data into a struct
% last modify: Lore 29/11
% Please before you had to run FastTool and do performance analysis 
% Turbine operations with Generate steady state operating conditions
load('SSData_IEA15MW')
SSData = struct();

SSData.U = WindSpeed;
SSData.PR = ElectricalPower;
SSData.T = RotorThrust;
SSData.Om = RotorSpeed;
SSData.B = PitchAngle;
SSData.TSR = TipSpeedRatio;

%% 1) Steady-state operating points
% Load input data
% load SSDataVS2PS1
% SSData = SSDataVS2PS1; clear SSDataVS2PS1

figure
subplot(221)
plot(SSData.U,SSData.PR/1e6,'o-','LineWidth',1.5); 
hold on, grid on
xlabel('Wind speed [m/s]'); ylabel('Aerodynamic power [MW]')
subplot(222)
plot(SSData.U,SSData.T/1e3,'o-','LineWidth',1.5); 
hold on, grid on
xlabel('Wind speed [m/s]'); ylabel('Aerodynamic thrust [kN]')
subplot(223)
plot(SSData.U,SSData.Om,'o-','LineWidth',1.5);
hold on, grid on
xlabel('Wind speed [m/s]'); ylabel('Rotor speed [rpm]')
subplot(224)
plot(SSData.U,SSData.B,'o-','LineWidth',1.5); 
hold on, grid on
xlabel('Wind speed [m/s]'); ylabel('Blade pitch [°]')

%% 2)
%load 'ASurf_IEA15MW'
load('Performance_ext.mat')
ASurf.B = Rotor_Pitch;
ASurf.Cp = Rotor_cP;
ASurf.Cq = Rotor_cQ;
ASurf.Ct = Rotor_cT;
ASurf.TSR = Rotor_Lamda';

Cp = ASurf.Cp; Cp(Cp<0) = NaN; %Blank coefficients where Cp<0
Ct = ASurf.Ct; Ct(ASurf.Cp<0) = NaN;

col = parula(10);

figure
subplot(1,2,1)
contourf(ASurf.TSR,ASurf.B,Cp);
hold on
plot(SSData.TSR,SSData.B,'x-r','linewidth',1.5)
colormap(col); caxis([0 0.5]); colorbar
xlabel('Tip-speed ratio [-]'); ylabel('Blade pitch [°]')
xlim([4 20]); ylim([-1 25])
subplot(1,2,2)
contourf(ASurf.TSR,ASurf.B,Ct);
hold on
plot(SSData.TSR,SSData.B,'x-r','linewidth',1.5)
colormap(col); caxis([0 1]); colorbar
xlabel('Tip-speed ratio [-]'); ylabel('Blade pitch [°]')
xlim([4 20]); ylim([-1 25])

% Additional data
rho = 1.225;
R   = 120;

%% ------------------------------------------------------------------------
% Torque coefficient from power coefficient (C_Q = C_P / \lambda)
for jj = 1:length(ASurf.B)
    ASurf.Cq(jj,:) = ASurf.Cp(jj,:)./ASurf.TSR;
end

% Torque from power (Q = P / RPM)
SSData.Q = SSData.PR./(SSData.Om*2*pi/60);

% Steady-state conditions
Cq0 = interp2(ASurf.TSR,ASurf.B,ASurf.Cq,SSData.TSR,SSData.B);
Q0  = 0.5*rho*pi*R^3*Cq0.*SSData.U.^2;
Ct0 = interp2(ASurf.TSR,ASurf.B,ASurf.Ct,SSData.TSR,SSData.B);
T0  = 0.5*rho*pi*R^2*Ct0.*SSData.U.^2;

%% ------------------------------------------------------------------------
% Sensitivities
% Gradient
dTSR = ASurf.TSR(2)-ASurf.TSR(1);
dB   = ASurf.B(2)-ASurf.B(1);

[pQTSR,pQB] = gradient(ASurf.Cq,dTSR,dB*pi/180);
[pTTSR,pTB] = gradient(ASurf.Ct,dTSR,dB*pi/180);

pQTSR0 = interp2(ASurf.TSR,ASurf.B,pQTSR,SSData.TSR,SSData.B);
pQB0   = interp2(ASurf.TSR,ASurf.B,pQB,SSData.TSR,SSData.B);

pTTSR0 = interp2(ASurf.TSR,ASurf.B,pTTSR,SSData.TSR,SSData.B);
pTB0   = interp2(ASurf.TSR,ASurf.B,pTB,SSData.TSR,SSData.B);

%-------------------------------------------
% SENSITIVITIES

% Rotor speed to rotor torque sensitivity
% K_{\omega Q} = \frac{Q_{r,0}}{\omega_0} \frac{\partial C_Q}{\partial\lambda}\biggl|_0 \frac{\l_0}{C_{Q,0}}
KOmQ = (Q0./(SSData.Om*2*pi/60)) .* (pQTSR0) .* (SSData.TSR./Cq0);

% Wind speed to rotor torque sensitivity
% K_{U Q} = \frac{Q_{r,0}}{U_0} \biggl(2 - \frac{\partial C_Q}{\partial\lambda}\biggl|_0 \frac{\l_0}{C_{Q,0}}\biggr)
KUQ = (Q0./(SSData.U)) .* (2 - (pQTSR0).*(SSData.TSR./Cq0));

% Rotor collective pitch angle to rotor torque sensitivity
% K_{\beta Q} = \frac{1}{2} \rho \pi R^3 U_0^2 \frac{\partial C_Q}{\partial\beta}\biggl|_0
KBQ = 0.5*rho*pi*R^3*SSData.U.^2.*pQB0;

% Rotor speed to rotor thrust sensitivity
% K_{\omega T} = \frac{T_0}{\omega_0} \frac{\partial C_T}{\partial\lambda}\biggl|_0 \frac{\l_0}{C_{T,0}}
KOmT = (T0./(SSData.Om*2*pi/60)).* (pTTSR0).* (SSData.TSR./Ct0);

% Wind speed to rotor thrust sensitivity
% K_{U T} = \frac{T_0}{U_0} \biggl(2 - \frac{\partial C_T}{\partial\lambda}\biggl|_0 \frac{\l_0}{C_{T,0}}\biggr)
KUT = (T0./(SSData.U)) .* (2 - (pTTSR0).*(SSData.TSR./Ct0));

% Rotor collective pitch angle to rotor thrust sensitivity
% K_{\beta T} = \frac{1}{2} \rho \pi R^2 U_0^2 \frac{\partial C_T}{\partial\beta}\biggl|_0
KBT = 0.5*rho*pi*R^2*SSData.U.^2.*pTB0;


%% ------------------------------------------------------------------------
% Figures
U = SSData.U;

%
figure
subplot(321)
plot(U,KOmQ,'-','LineWidth',3)
hold on
title('Torque')
ylabel('K_{\omega{Q}} [Nm/rad/s]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
subplot(323)
plot(U,KUQ,'-','LineWidth',3)
hold on
ylabel('K_{UQ} [Nm/m/s]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
subplot(325)
plot(U,KBQ,'-','LineWidth',3)
hold on
ylabel('K_{\beta{Q}} [Nm/rad]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
subplot(322)
plot(U,KOmT,'-','LineWidth',3)
hold on
title('Thrust')
ylabel('K_{\omega{T}} [N/rad/s]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
subplot(324)
plot(U,KUT,'-','LineWidth',3)
hold on
ylabel('K_{UT} [N/m/s]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
subplot(326)
plot(U,KBT,'-','LineWidth',3)
hold on
ylabel('K_{\beta{T}} [N/rad]')
xlabel('Wind speed [m/s]')
set(gca,'FontSize',9)
xlim([3 25])
