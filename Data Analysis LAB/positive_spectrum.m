function [pos_spec,f_pos]=positive_spectrum(neg_spec,fs)

% The function takes the spectrum of a signal defined also in negative
% frequencies and gives back the spectrum defined only in positive
% frequencies.
%
%   [pos_spec,f_pos]=positive_spectrum(neg_spec,fs)
%
% Inputs:
%   - neg_spec: spectra of the signal defined both in positive and 
%               negative frequencies (spectra computed with fft, for ex). 
%               It's a matrix whose columns are the spectra
%   - fs: sampling frequency (optional, needed if you want f_pos)
%
% Outputs:
%   - pos_spec: spectra of the signal defined only in positive frequencies
%               (again, it's a matrix whose columns are the spectra)
%   - f_pos: positive frequencies over which the new spectrum is defined;
%            requires fs as input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% informations about the inputs

N=length(neg_spec);

if nargin>1
    
    % creating frequency vector

    df=fs/N;
    if floor(N/2)==N/2  % even number of frequencies
        f_pos=0:df:(N/2*df);
    else    % odd number of frequencies
        f_pos=0:df:((N-1)/2)*df;
    end

end

% computing normalised spectra

if floor(N/2)==N/2  % even number of frequencies
    pos_spec(1,:)=neg_spec(1,:);
    pos_spec(2:N/2,:)= neg_spec(2:N/2,:).*2;
    pos_spec(N/2+1,:)=neg_spec(N/2+1,:);
else    % odd number of frequencies
    pos_spec(1,:)=neg_spec(1,:);
    pos_spec(2:(N+1)/2,:)= neg_spec(2:(N+1)/2,:).*2;
end

end
