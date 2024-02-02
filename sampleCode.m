%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample codes highlighting the usage of provided functions
% from the Doubly-Dispersive Channel Communications Package 
% - Eric Hyeon-Seok Rou, 2024.
% [https://arxiv.org/abs/2401.07700] 
close all; clear all; clc;

%% Simulation parameters
carrierFreq = 4e9;  % Carrier frequency in [Hz]
carrierWavelength = physconst("lightspeed")/carrierFreq;
bandwidth = 1e6;    % Bandwidth in [Hz]
samplingFreq = 2*bandwidth;
%
bps = 4;    % Bits per symbol
M = 2^bps;  % Symbol constellation cardinality
N = 2^8;    % Number of subcarriers
NT = 1;     % Number of symbol periods

%% Define channel statistics
maxDelay = 16;       % Maximum delay index
maxDoppler = 2;      % Maximum normalised digital velocity
maxRange = (maxDelay/samplingFreq)/2*physconst("lightspeed"); % Maximum unambiguous range in [m] 
maxVel = (maxDoppler/N*samplingFreq)/2*carrierWavelength;     % Maximum unambiguous velocity in [m/s]

%% Generate Tx Symbols
X_i = randi([0 M-1], N, NT);                   % Index of transmit symbols
X = qammod(X_i, M, "UnitAveragePower", true);  % Complex transmit symbols

%% Modulate with AFDM 
guardWidth = 4;     % Guard width of the AFDM waveform
if (2*(maxDoppler + guardWidth)*(maxDelay + 1)) + (maxDelay) > N
    fprintf("Integer orthogonality of the AFDM is not satisfied!\n");
end
c1_AFDM = (2*(maxDoppler + guardWidth) + 1)/(2*N); % Satisfying the orthogonality condition
c2_AFDM = 1/(N^2*pi);                              % Sufficiently small irrational number
%
S_AFDM = AFDMmod(X, c1_AFDM, c2_AFDM);             % Modulate symbol X with AFDM

%% Modulate with OTFS
% Arbitrarily set up the K x L delay-Doppler grid size, such that K*L = N 
facN = factor(N); facIn = min(max(round((length(facN)/2) + randn(1)), 1), length(facN)); 
K = prod(facN(1:facIn));
L = prod(facN(facIn+1:end));
X_OTFS = reshape(X, K, L, []);
%
S_OTFS = OTFSmod(X_OTFS);                          % Modulate symbol X with OTFS

%% Demodulation after AWGN channel
SNRdB = 30;
%
X_est_AFDM = AFDMdemod(awgn(S_AFDM, SNRdB), c1_AFDM, c2_AFDM);
X_est_OTFS = OTFSdemod(awgn(S_OTFS, SNRdB), [K, L]);


%% Generate a Doubly-Dispersive Channel Realization with Random Paths
P = 3;                                            % Number of scattering paths 
%
pathRange = maxRange*rand(1, P);                  % random scatter range in [m]
pathVel = maxVel*cos((2*randn(1,P) - 1)*pi);      % random scatter relative velocity in [m/s]
%
pathDelays = 2*pathRange/physconst("lightspeed"); % path delay in [s]
pathDopplers = 2*pathVel/carrierWavelength;       % path Doppler shift in [Hz]
%
fadingCoefs = (randn(1,P)+1j*randn(1,P))/sqrt(2); % channel fading coefficients (Rayleigh)
% circular convolutional channel matrix with non-chirp cyclic prefix (just CP)
[H_noCCP] = cconvLTVChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, false);
% imagesc(abs(H_noCCP));
%
% circular convolutional channel matrix with chirp-cyclic prefix 
[H_wCCP] = cconvLTVChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, true, c1_AFDM);
imagesc(abs(H_wCCP));

R_AFDM = H_wCCP  * S_AFDM;
R_OTFS = H_noCCP * S_OTFS;
