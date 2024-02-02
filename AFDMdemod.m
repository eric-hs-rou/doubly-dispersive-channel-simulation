function [X] = AFDMdemod(S, c1, c2)
% X = AFDMdemod(S, c1, c2) AFDM demodulation with defined chirp frequencies
%   Inputs:
%      - S  : N x T sampled sequence of a AFDM modulated signal
%      - c1 : central digital frequency of first AFDM chirp
%      - c2 : central digital frequency of second AFDM chirp
%   Output:
%      - X : N x T matrix of AFDM demodulated samples

sizeS = size(S);
N = sizeS(1);

L1 = diag(exp(-1j*2*pi*c1*((0:N-1).^2)));
L2 = diag(exp(-1j*2*pi*c2*((0:N-1).^2)));
F = dftmtx(N)/sqrt(N);

A = L2*F*L1;

X = A*S;

end