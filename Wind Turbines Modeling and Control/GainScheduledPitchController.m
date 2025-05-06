% Gain scheduling and certification
% Run FastTool and do performance analysis 
% Turbine operations with Generate steady state operating conditions

%% Gain Scheduling
Vr = 10.37; % Rated wind speed TO BE CHECKED
indice_vr = find(U >= Vr); % Indice che serve per andare a prendere i valori ottimali delle sensitivity function
%indice_vr = indice_vr(1);
indice_vr = find(U >= 11);


% wind_start = 10;
step_U = U(2)-U(1);
% index = wind_start/step_U+1:height(U);

figure("Name", "Rotor collective pitch angle to rotor torque sensitivity")  % "Pitch Abgle to Wind Speed")
yyaxis left
plot(U(indice_vr),KBQ(indice_vr),'o-','LineWidth',1)
ylabel('K_{\beta{Q}} [Nm/rad]')
xlabel('Wind speed [m/s]')

% figure
% plot(dCp_dPitch(:,NewIdxTSR)/Qg§(NewIdxTSR))

hold on
yyaxis right
plot(U(indice_vr),PitchAngle(indice_vr),'o-', 'LineWidth',1);
xlabel('Wind speed [m/s]'); ylabel('Blade pitch [deg]');

%% Plant Model - effect of the steady-state operating conditions
P = polyfit(PitchAngle(indice_vr), KBQ(indice_vr), 1);
KBQ_fit = polyval(P, PitchAngle(indice_vr));
%KBQ_fit = P(1) * PitchAngle(index) + P(2);

figure
plot(PitchAngle(indice_vr), KBQ(indice_vr),'b-o', PitchAngle(indice_vr),KBQ_fit,'r-','LineWidth',1)
ylabel('K_{\beta{Q}} [Nm/rad]')
xlabel('Blade pitch [deg]')
legend('Real Data', 'Linear fitting', 'Location','best')


% Guardo Gain scheduling E06 per le notazioni
KOmQ_0 = KOmQ(indice_vr); 
KBQ_0 = KBQ(indice_vr); 

BladePitch = PitchAngle(indice_vr)*pi/180;
P1 = polyfit(BladePitch - BladePitch(1), KBQ(indice_vr),1);
C1 = P1(1)/P1(2)
C2 = 1;


%% Gain scheduling depending on wind speed designed
% from U = 10.5 to U = 11.2
set_U = [10.5:0.1:11.2,12];
Kpd = [-3.3338 -3.2206 -2.2958 -2.2352 -1.7765 -1.7595 -1.7427 -1.5796 -1.0415];
Kid = [-0.8105 -0.7795 -0.5514 -0.5362 -0.4281 -0.4221 -0.4161 -0.3831 -0.2618];

for ii = 1:1:length(Kpd)
wind_speed_des = set_U(ii); % scelto da noi e che possiamo cambiare
kp_des = Kpd(ii);
ki_des = Kid(ii);
Beta_i = PitchAngle(int64(wind_speed_des/step_U));
disp(['Angolo Beta_i : ', num2str(Beta_i), ' per U = ', num2str(wind_speed_des)]);
% calcolo la funzione eta(beta) scegliendo come riferimento 11 m/s

eta = @(Beta) (C1*(Beta_i*pi/180)+1)./(C1*(Beta*pi/180)+1);

% Set up for plots
Beta_linspace = linspace(0,25,2500);
legenda = ['Designed wind speed: ', num2str(wind_speed_des)];
figure("Name",'Gain Scheduling')
subplot(2,1,1)
plot(Beta_linspace, eta(Beta_linspace), 'LineWidth',1.5)
hold on
xline(Beta_i, 'k--')
xlabel('Blade Pitch [°]');
ylabel('Gain scheduling (\eta(\beta)) [-]');
legend(legenda)
hold off

subplot(2,1,2)
plot(Beta_linspace, eta(Beta_linspace) * Kpd(ii),'r' ,'LineWidth',1.5)
hold on
plot(Beta_linspace, eta(Beta_linspace) * Kid(ii),'b', 'LineWidth',1.5)
xline(Beta_i, 'k--')
xlabel('Blade Pitch [°]');
ylabel('Gain Scheduling [\delta\beta]');
legend('K_p','K_i','\beta_i','Location','best')
end

