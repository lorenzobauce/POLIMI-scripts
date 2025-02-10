function coh=coherence(S_xx,S_yy,S_xy,option)

% Compute the coherence function between two signals x and y.
%
%   coh=coherence(S_xx,S_yy,S_xy,option)
%
% Inputs:
%   - S_xx: power spectrum of the first signal (also G_xx is good)
%   - S_yy: power spectrum of the second signal (also G_yy is good)
%   - S_xy: cross spectrum (also G_xy is good)
%   - option: "1" if you are giving G_xx, G_yy, G_xy (spectra defined 
%             only in positive frequencies) as inputs, instead of
%             S_xx, S_yy and S_xy (default: "0")
%
% Output:
%   - coh: coherence function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    % convert spectra only in positive frequencies
    G_xx=positive_spectrum(S_xx);
    G_yy=positive_spectrum(S_yy);
    G_xy=positive_spectrum(S_xy);
else
    if option~=1
        % convert spectra only in positive frequencies
        G_xx=positive_spectrum(S_xx);
        G_yy=positive_spectrum(S_yy);
        G_xy=positive_spectrum(S_xy);
    else
        G_xx=S_xx;
        G_yy=S_yy;
        G_xy=S_xy;
    end
end

% compute coherence

coh=abs(G_xy).^2./(G_xx.*G_yy);

end