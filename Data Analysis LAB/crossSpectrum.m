function [Sxy,f]=crossSpectrum(x,y,T,fs,ovrl)

% The function computes the cross spectrum of two signals x and y, 
% averaging many cross spectra of subrecords of the signals themselves. 
% Hanning windowing is used. To obtain the power spectrum, just input the 
% same signal twice.
%
%   [Sxy,f]=crossSpectrum(x,y,T,fs,ovrl)
%
% Inputs:
%   - x: vector containing the first sampled signal
%   - y: vector containing the second sampled signal
%   - T: time width of the subrecords we want to use for averaging (notice
%        that this determines the frequency resolution of the cross 
%        spectrum: df=1/T)
%   - fs: sampling frequency of the signal records
%   - ovrl: (optional) overlap of the subrecords (scalar from 0 to <1, default: 0)
%
% Outputs:
%   - f: frequencies over which the cross spectrum is computed, as a vector
%   - Sxy: cross spectrum, as a vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check

if nargin>3
    if ovrl<0 || ovrl>=1
        error('ERROR in subrec.m: the overlap must be from 0 to <1')
    end
else
    ovrl=0;
end

if length(x)~=length(y)
    error('ERROR in crossSpectrum.m: data vectors must have the same length')
end

% compute relevant parameters

dt=1/fs; % sampling time
df=1/T; % frequency resolution
N_sub=ceil(T/dt);    % number of samples in each subrecord

% compute the frequency axis over which the spectrum is defined

f=-df*N_sub/2:df:df*(N_sub/2-1);

% extract subrecords from the signal 

[sub_mat_x]=subrec(x,T,fs,ovrl);
[sub_mat_y]=subrec(y,T,fs,ovrl);

% window each subrecord using Hanning windowing, to reduce leakage

w=hann(N_sub,'periodic');
sub_mat_x=sub_mat_x.*w;
sub_mat_y=sub_mat_y.*w;

% now compute the spectrum AND the cross spectrum for each windowed
% subrecord and average them

X_mat_x=fft(sub_mat_x); % X_mat is a matrix whose columns are the spectra of the corresponding subrecords
X_mat_y=fft(sub_mat_y);
X_mat_x=X_mat_x./N_sub; % rescaling the spectrum (because the matlab funtion fft is defined like shit)
X_mat_y=X_mat_y./N_sub;
Sxy_mat=conj(X_mat_x).*(X_mat_y); % matrix containing the cross spectra as columns
Sxy=mean(Sxy_mat,2);  % averaged cross spectrum

end