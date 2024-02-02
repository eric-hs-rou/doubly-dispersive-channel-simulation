function [H, H_p] = cconvLTVChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, includeCCP, varargin)
% H = circConvDDChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, modType, [c1])
% Generates the circular convolutional matrix form of the doubly-dispersive channel, given scattering point parameters
%   Inputs:
%      - N            : size of channel matrix
%      - pathDelays   : P x 1 vector of path delays
%      - pathDopplers : P x 1 vector of path Doppler shifts
%      - fadingCoefs  : P x 1 vector of channel fading coefficients
%      - samplingFreq : sampling frequency of the system
%      - includeCCP   : include chirp-cyclic prefix effect (true/false)
%      - c1           : first AFDM chirp frequency (only required if "CCP" is "true")
%   Output:
%      - H   : N x N Full circular convolutional channel matrix
%      - H_p : N x N x P circular convolutional channel matrix per path

if length(pathDelays) == length(pathDopplers) && length(pathDopplers) == length(fadingCoefs)
    P = length(pathDelays);
else
    error("Invalid number of target parameters");
end

normDelays = round(pathDelays*samplingFreq);   % normalised integer delay index
normDopplers = N*pathDopplers/samplingFreq;    % normalised digital Doppler shift

W = diag(exp(-1j*2*pi*(0:N-1)/N));
Pi = [zeros(1, N-1), 1; eye(N-1), zeros(N-1, 1)];

if includeCCP
    ccpf = varargin{1};
else
    ccpf = 0;
end

H_p = zeros(N, N, P);
for p = 1:P
    CPP_p = diag( [exp(-1j*2*pi*ccpf*((N^2)-2*N*( normDelays(p) - (0:(normDelays(p)-1) )))), ones(1, N - normDelays(p)) ] ) ;
    H_p(:,:,p) = fadingCoefs(p) * CPP_p * W^normDopplers(p) * Pi^normDelays(p);
end
H = sum(H_p, 3);


end