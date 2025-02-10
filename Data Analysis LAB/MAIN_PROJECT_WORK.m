clear
close all
clc

%% PLOTS CONTROL

plot_time=1;        % signals in time domain
plot_freq=1;        % signals in frequency domain
plot_powcross=1;    % power and cross spectra
plot_coh=1;         % coherence function
plot_tf=1;          % transfer function estimation (H1 and H2)
plot_modrec=1;      % reconstruction of estimated transf function

%% TRANSFER FUNCTION ESTIMATION FROM EXPERIMENTAL DATA

%% Load experimental data

% load data: we expect 2 records: random input (white noise) (first in the 
% data file) and its output (second in the data file); we also need to know
% the sampling frequency and the sensitivity.

[data_file,path]=uigetfile;
cd(path);
load(data_file);

%% Set/compute test parameters

% change names to the data so that the code runs

x=Dati(:,2);    % input is called x
y=Dati(:,1);    % output is called y
fsamp=2048;    % sampling frequency
sens_in=10*1e-3;      % input sensitivity [V/A]
sens_out=102*1e-3/9.81;     % output sensitivity [V/g]

% compute relevant parameters

n=max(size(y));    % number of samples
dt=1/fsamp;        % sampling time
t_end=(n-1)*dt;    % final time istant
t=0:dt:t_end;      % time vector

%% Data conversion

% convert from voltage to force and acceleration using sensitivities

x=x./sens_in;   % now x is in [A]
y=y./sens_out;  % now y is in [m/s^2]


%% Plot signals in time domain

if plot_time==1

figure

subplot(2,1,1)  % plot input force
plot(t,x)
grid on
title('Input force in time domain')
xlabel('t [s]')
ylabel('x [A]')
axis tight

subplot(2,1,2)  % plot output
plot(t,y)
grid on
title('Output in time domain')
xlabel('t [s]')
ylabel('y [m/s^2]')
axis tight

end

%% Signals in frequency domain

% compute the spectra of the signals

X=fft(x);
Y=fft(y);

% convert in only positive frequencies

[X_pos,f_pos]=positive_spectrum(X,fsamp);
Y_pos=positive_spectrum(Y,fsamp);

% plot input force in frequency domain

if plot_freq==1

figure

subplot(2,1,1)  % amplitude
semilogy(f_pos,abs(X_pos))
grid on
title('Input force in frequency domain')
xlabel('f [Hz]')
ylabel('|X|')
axis tight

subplot(2,1,2)  % phase
plot(f_pos,rad2deg(angle(X_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠X  [deg]')
axis tight

% plot output

figure

subplot(2,1,1)  % amplitude
semilogy(f_pos,abs(Y_pos))
grid on
title('Output in frequency domain')
xlabel('f [Hz]')
ylabel('|Y|')
axis tight

subplot(2,1,2)  % phase
plot(f_pos,rad2deg(angle(Y_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠Y  [deg]')
axis tight

end

%% Averaging parameters definition

% overlap and time duration of the subrecords used to average the power and
% cross spectra

T=10;       % assign time width of the subrecords
overlap=0.66;  % assign overlap of the subrecords

%% Power spectra and cross spectra (actually PSD)

% compute power and cross spectra

[Sxx,f_T]=crossSpectrum(x,x,T,fsamp,overlap);
Syy=crossSpectrum(y,y,T,fsamp,overlap); 
Sxy=crossSpectrum(x,y,T,fsamp,overlap);

% normalise for delta_f to obtain PSD: signals are random, so we need to
% take in account for leakage

delta_f=1/T;
Sxx=Sxx./delta_f;
Syy=Syy./delta_f;
Sxy=Sxy./delta_f;

% convert for only positive frequencies

[Gxx,f_T_pos]=positive_spectrum(Sxx,fsamp);
Gyy=positive_spectrum(Syy);
Gxy=positive_spectrum(Sxy);

% plot them

if plot_powcross==1

figure   % PSD of input
semilogy(f_T_pos,Gxx) 
xlabel('f [Hz]')
ylabel('G_x_x [A^2/Hz]')
grid on
title('PSD of the input')
axis tight

figure  % PSD of output
semilogy(f_T_pos,Gyy)
xlabel('f [Hz]')
ylabel('G_y_y [(m/s^2)^2/Hz]')
grid on
title('PSD of the output')
axis tight

figure  % cross spectrum

subplot(2,1,1)  % magnitude
semilogy(f_T_pos,abs(Gxy))
xlabel('f [Hz]')
ylabel('|G_x_y|')
grid on
title('Cross spectrum density')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(Gxy)))
xlabel('f [Hz]')
ylabel('∠G_x_y  [deg]')
grid on
axis tight

end

%% Coherence function

% compute coherence function

coh=coherence(Sxx,Syy,Sxy);

% plot it

if plot_coh==1

figure
plot(f_T_pos,coh);
grid on
xlabel('f [Hz]')
ylabel('\gamma^2_x_y')
title('Coherence function')
axis tight

end

%% Transfer function estimation

% compute H1 and H2

Gyx=conj(Gxy);
H1_pos=Gxy./Gxx;
H2_pos=Gyy./Gyx;

% plot H1 and H2

if plot_tf==1

figure  % H1

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H1_pos))
grid on
title('H_1 estimator')
xlabel('f [Hz]')
ylabel('|H_1|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_1  [deg]')
axis tight

figure  % H2

subplot(2,1,1)  % amplitude
semilogy(f_T_pos,abs(H2_pos))
grid on
title('H_2 estimator')
xlabel('f [Hz]')
ylabel('|H_2|')
axis tight

subplot(2,1,2)  % phase
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('∠H_2  [deg]')
axis tight

end

%% Comparisons - H1, H2, and coherence

% plot H1, H2 and the coherence together

if plot_tf==1

figure

subplot(3,1,1)  % coherence
plot(f_T_pos,coh,'k')
grid on
title('Coherence')
xlabel('f [Hz]')
ylabel('\gamma_x_y^2')
axis tight

subplot(3,1,2)  % H1 and H2 amplitude
semilogy(f_T_pos,abs(H1_pos))
hold on
semilogy(f_T_pos,abs(H2_pos))
grid on
title('Comparison between H_1 and H_2')
xlabel('f [Hz]')
ylabel('amplitude')
legend('H_1','H_2')
axis tight
hold off
    
subplot(3,1,3)  % phase
plot(f_T_pos,rad2deg(angle(H1_pos)))
hold on
plot(f_T_pos,rad2deg(angle(H2_pos)))
grid on
xlabel('f [Hz]')
ylabel('phase [deg]')
legend('H_1','H_2')
axis tight
hold off

end

%% Extraction of "relevant" part of H1 and H2

% plot coherence along with H1 and H2 in the same plot

figure

yyaxis left
ylabel('amplitude')
semilogy(f_T_pos,abs(H1_pos),'b-')
hold on
semilogy(f_T_pos,abs(H2_pos),'r-')

yyaxis right
ylabel('\gamma_x_y^2')
plot(f_T_pos,coh,'k')

grid on
title('Choose the freqneucy range to analyse')
xlabel('f [Hz]')
legend('H_1','H_2','coherence')
axis tight
hold off

% choose frequency range we want to consider (where coherence is high)

[select,~]=ginput(2);
f_initial=select(1);                                  % initial frequecy value
f_final=select(2);                                    % final frequency value

f_initial_indx=find(f_T_pos>=f_initial,1,'first');    % initial frequecy index
f_final_indx=find(f_T_pos>=f_final,1,'first');        % final frequecy index

% redefine everything inside this frequency range

f_range=f_T_pos(f_initial_indx:f_final_indx);
H1_range=H1_pos(f_initial_indx:f_final_indx);
H2_range=H2_pos(f_initial_indx:f_final_indx);
coh_range=coh(f_initial_indx:f_final_indx);

% plot to check

figure

yyaxis left
semilogy(f_range,abs(H1_range),'b-','LineWidth',1.5)
hold on
semilogy(f_range,abs(H2_range),'r-','LineWidth',1.5)
ylabel('amplitude')

yyaxis right
plot(f_range,coh_range,'k--')
ylabel('\gamma_x_y^2')

grid on
title('H_1, H_2 and coherence in the selected range')
xlabel('f [Hz]')
legend('H_1','H_2','coherence')
axis tight
hold off

%% MODAL RECONSTRUCTION

% We want to reconstruct modal parameters (in particular natural 
% frequencies and damping) of our system by looking at experimental data 
% (a matrix of experimental transfer functions). This matrix is a A x R
% matrix, where A is the length of the frequency vector, while R is the
% number of tranfer functions (each r-th column is associated to a tranfer 
% function). In our case R=1: we only have a vector containing H1.

% domains definition and change variables names (to make the two codes
% match)

f=f_range;      % frequencies [Hz]
w=2*pi*f;       % pulsations [rad/s]
G_exp=H1_range; % our experimental data is H1, in the range we provided
n_in=1;         % we have 1 input
n_out=1;        % we have 1 output

%% Finding natural frequencies

% adjust the prominence for findpeaks (increase it if it finds too many
% peaks)

prominence=1;

% find natural frequencies (as peaks of H1)

[~,indx]=findpeaks(log(abs(G_exp(:,1))),'MinPeakProminence',prominence);
f_0=f(indx);    % natural frequencies

% visual check (adjust if necessar)

figure

yyaxis left
semilogy(f_range,abs(H1_range),'b-','LineWidth',1.5)
hold on
semilogy(f_range,abs(H2_range),'r-','LineWidth',1.5)
semilogy(f_0,abs(G_exp(indx,1)),'v','MarkerFaceColor','b','MarkerSize',10)
ylabel('amplitude')

yyaxis right
plot(f_range,coh_range,'k--')
ylabel('\gamma_x_y^2')

grid on
title('Natural frequencies in the selected range')
xlabel('f [Hz]')
legend('H_1','H_2','natural frequencies','coherence')
axis tight
hold off

fprintf(['Natural frequencies found in the experimental data [Hz]:\n', num2str(f_0,5),'\nCorresponding to indices:\n', num2str(indx'),'\n\n'])

% allocate the informations we found

n_modes=length(indx);   % how many modes
w_0=w(indx);            % natural pulsations

%% Define optimization parameters

% "width" of our approximation (in terms of n° of indices)

opt_length=inputdlg(['Choose optimisation range, in terms of range of indices (consider that first natural frequency is at at index ', num2str(indx(1)),')']);
opt_length=str2double(opt_length);

%% Modal reconstruction/optimisation procedure

% we perform the optimization for each mode (resonance peak)

f_opt_mat=zeros(opt_length+1,n_modes);           % first output of the for cycle: frequency vectors (as columns) where we perform the optimization (1 for each mode)
modal_param_mat=zeros(2+3*n_in*n_out,n_modes);   % second output of the for cycle: matrix that contains (each column) the extimated X (containing modal
                                                 % parameters w_k,xi_k, A_k, RL_k and RH_k); 1 column for each mode
for kk=1:n_modes

    % compute the frequency range around the k-th mode where we perform the
    % optimization and allocate it in the correponding column of f_opt_mat
    
    indx_in=indx(kk)-ceil(opt_length/2);
    indx_end=indx(kk)+floor(opt_length/2);
    f_opt=f(indx_in:indx_end);

    f_opt_mat(:,kk)=f_opt';
    
    % Isolate the experimental transfer funct matrix around the k-th peak.
    % G_exp_k contains U rows and R columns (each column contains a transf
    % funct, defined over the "restricted" U frequencies)
    
    U=length(f_opt);
    R=min(size(G_exp));
    
    G_exp_k=G_exp(indx_in:indx_end,:);
    
    % define the initial guessess vector for the unknown parameters
    
    w_k_guess=w_0(kk);   % k-th natural pulsation (take it close to the one we know from data)
    
    phase_der=( angle(G_exp(indx(kk),1))-angle(G_exp(indx(kk)-1,1)) ) / ( 2*pi*f(indx(kk))-2*pi*f(indx(kk)-1) ); % adimensional damping estimated with phase derivative method
    xi_k_guess=-1/(w_k_guess*phase_der);
    
    A_k_guess=zeros(1,R);      % vector of coefficients A_k (each element refers to the r-th trans funct):
                               % each element is choosen as minus the imaginary part of the experimental transfer
                               % transf funct evaluated at the k-th natural freq, multiplied by xi_k and w_k
    for rr=1:R
        [~,indx_resonance]=findpeaks(abs(G_exp_k(:,rr)),'NPeaks',1);
        A_k_guess(rr)=-imag(G_exp_k(indx_resonance,rr))*2*xi_k_guess*w_k_guess^2;
    end
    
    RL_k_guess=zeros(1,R)+eps*1j*ones(1,R); % left-residual coefficients
    RH_k_guess=zeros(1,R)+eps*1j*ones(1,R); % right-residual coefficients
    
    X_guess=[w_k_guess,xi_k_guess,A_k_guess,RL_k_guess,RH_k_guess]';   % initial guess parameters vector
    
    % compute the error function that only depends on parameters X to be found
    % by lsqnonlin (see parameterfun for more information)
    
    err=@(X) parameterfun(X,f_opt,G_exp_k);
    
    % minimise with matlab function
    
    options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
    X=lsqnonlin(err,X_guess,[],[],options);
    
    modal_param_mat(:,kk)=X;
end
    
%% Checking the extimated parameters
    
if plot_modrec==1

% plot amplitudes

figure
mi=1;

for mm=1:n_in   
    for ii=1:n_out
        
        % plot experimental data

        subplot(n_in,n_out,mi)
        semilogy(f,abs(G_exp(:,mi)),'LineWidth',1.5)
        hold on

        % plot numerical, reconstructed transfer functions
        
        for kk=1:n_modes
            G_num=@(Omega) (modal_param_mat(mi+2,kk)./(-Omega.^2+1j*2*modal_param_mat(2,kk)*modal_param_mat(1,kk)*Omega+modal_param_mat(1,kk)^2)) + (modal_param_mat(mi+R+2,kk)./(Omega.^2)) + (modal_param_mat(mi+2*R+2,kk));
            ampl=abs(G_num(2*pi*f_opt_mat(:,kk)));
            semilogy(f_opt_mat(:,kk),ampl,'LineWidth',1,'Color','r','Marker','.')           
        end

        % plot grid, axis, etc

        grid on
        xlabel('f [Hz]')
        ylabel('amplitude')        
        mi=mi+1;
        title('Magnitude modal reconstruction')
        axis tight
        hold off

    end
end

% plot phases

figure
mi=1;

for mm=1:n_in 
    for ii=1:n_out

        % plot experimental data

        subplot(n_in,n_out,mi)
        plot(f,angle(G_exp(:,mi)), 'LineWidth',1.5)
        hold on

        % plot numerical, reconstructed transfer functions
        
        for kk=1:n_modes
            G_num=@(Omega) (modal_param_mat(mi+2,kk)./(-Omega.^2+1j*2*modal_param_mat(2,kk)*modal_param_mat(1,kk)*Omega+modal_param_mat(1,kk)^2)) + (modal_param_mat(mi+R+2,kk)./(Omega.^2)) + (modal_param_mat(mi+2*R+2,kk));
            phase=angle(G_num(2*pi*f_opt_mat(:,kk)));
            plot(f_opt_mat(:,kk),phase,'LineWidth',1,'Color','r','Marker','.')           
        end

        % plot grid, axis, etc

        grid on
        xlabel('f [Hz]')
        yticks([-pi -pi/2 0 pi/2 pi])
        yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
        ylabel('phase [rad]')
        title('Phase modal reconstruction')
        axis tight
        hold off
        mi=mi+1;

    end
end

end

% print and save the estimated modal parameters

nat_freq_estim=modal_param_mat(1,:)./2/pi;  % natural frequencies
damp_estim=modal_param_mat(2,:); % adimensional damping coefficients

fprintf(['\nReconstructed natural frequencies [Hz]:\n', num2str(nat_freq_estim,5),'\n\n'])
fprintf(['Reconstructed adimensional dampings [-]:\n', num2str(damp_estim,5),'\n\n'])

%% SALVATAGGIO DATI
 
% sampling_frequencies=fsamp;
% natural_frequencies=f_0;
% natural_pulsations=w_0;
% estimated_dampings=damp_estim;
% H1=H1_range;
% indici_picchi=indx;
% input_sensitivity=sens_in;
% output_sensitivity=sens_out;
% unita_di_misura='input: ampere [A], output: accelerazione [m/s^2]';
% 
% %%
% save('circuit_mainVariables_VISE.mat','natural_frequencies','sampling_frequencies','natural_pulsations','estimated_dampings','H1','indici_picchi','input_sensitivity','output_sensitivity','unita_di_misura')
% save('circuit_allWorkspace_VISE.mat')

%%

%load("tutto.mat")