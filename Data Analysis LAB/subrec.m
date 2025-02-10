function [sub_mat]=subrec(data,T,fs,ovrl)

% The function takes a sampled signal and extracts subrecords of desired 
% duration from it, arranging them in the columns of a matrix. The 
% subrecords can also have an overlap.
%
%   [sub_mat]=subrec(data,T,fs,ovrl)
%
% Inputs:
%   - data: vector containing the sampleed signal
%   - T: time width of the subrecords we want to extract
%   - fs: sampling frequency of the sampled signal
%   - ovrl: (optional) overlap of the subrecords (scalar from 0 to <1, default: 0)
%
% Output:
%   - sub_mat: matrix containing the extracted subrecords as columns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check

if nargin>3
    if ovrl<0 || ovrl>=1
        error('ERROR in subrec.m: the overlap must be from 0 to <1')
    end
else
    ovrl=0;
end

% extracting relevant parameters

N=length(data); % number of samples in the signal
dt=1/fs;    % sampling time
N_sub=ceil(T/dt);    % number of samples in each subrecord
N_ovrl=round(ovrl*N_sub);  % overlapped samples
N_incr=N_sub-N_ovrl;    % increment of the "extracting" window
M=ceil((N-N_sub+1)/N_incr);   % number of expected subrecords

% extract all the subrecords and group them in a matrix:

sub_mat=zeros(N_sub,M);
for ii=1:M 
    % define first and last indices of the ii-th subrecord
    first_indx=(ii-1)*N_incr+1;
    fin_indx=first_indx+N_sub-1;
    
    % extract from the singnal and put in the matrix
    sub_mat(:,ii)=data(first_indx:fin_indx);
end

end