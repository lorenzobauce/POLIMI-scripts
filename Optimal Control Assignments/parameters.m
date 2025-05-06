%Script con le variabili globali usati nelle funzioni che danno le matrici
%dei sistemi lineari e non

%% VALORI ATTRITO ED ELASTICI STIMATI SPERIMENTALMENTE
global Jtheta_corr Kwire Cenc
Jtheta_corr=1.63e-5; % Fattore correttivo dell'inerzia rispetto l'asse z (theta)
Kwire=1.167e-4; % Ottenuto ottmizzando
Cwire=9.68e-5; % Ottenuto ottimizzando
Cenc=1.65e-5; % Sperimentalmente era stato trovato 7.5e-06, ma non era troppo attendibile; % to be checked

%% PARAMETRI ELETTRICI MOTORE
global K_tau R_tau Ctheta
K_tau=0.042;  % Motor back-emf constant
R_tau=8.4; % [Ohm] motor resistance circuit
C_tau=(K_tau^2)/R_tau; % Resistenza viscosa equivalente della veck-emf
Ctheta=C_tau+Cwire; % Corrisponde al 3e-4 trovato sperimentalmente

%% PARAMETRI DI INERZIA
global g Jrod Jenc Jtau mp r Lp

g=9.80665; % Accelerazione gravità (m/s^2)

% MOTORE E ATTACCO DEL MOZZO
Jtau=(0.40+0.06)*10^-5; %momento di inerzia da produttore (kg m^2)

% CORPO ENCODER come cilindro che ruota attorno ad un asse perpendicolare traslato dal baricentro
me=82.9*10^-3; %massa (kg)
Le=51.5*10^-3; %lunghezza (m)
Rin=3.15*10^-3; %raggio interno (m)
Rex=14.4*10^-3; %raggio esterno (m)
De=7.05*10^-3; %distanza asse rotazione dal baricentro (m)
Jenc=1/12*me*(3*(Rin^2+Rex^2)+Le^2)+me*(De^2);  %momento di inerzia (kg m^2)

% ASTA ROD come asta che ruota attorno ad un asse asse passante per un'estremità
mr=9*10^-3; %massa (kg)
Lrhat=101.6*10^-3; %lunghezza (m)
Dr=29.8*10^-3; %distanza asse rotazione dal baricentro (m)
Jrod=1/12*mr*Lrhat^2+mr*(Dr^2); %momento di inerzia (kg m^2)

% ASTA PENDOLO
mp=24.10*10^-3; %massa (kg)
Lp=128.8*10^-3;; %lunghezza (m)
r=85*10^-3;; %distanza dall'asse z del piano in cui oscilla (m)