function [X] = OTFSdemod(varargin)
% X = OTFSdemod(S, [K, L], pf) OTFS modulation
%   Inputs:
%      - S      : KL * T sample matrix
%      - [K, L] : Size of the OTFS delay-Doppler grid
%      - pf     : K x 1 receive pulse-shaping filter (rectangular pulse if not provided)
%   Output:
%      - X  : K x L x T tensor of OTFS demodulated symbols

if nargin == 2

    S = varargin{1};
    KL = varargin{2};
    K = KL(1);
    L = KL(2);

    sizeS = size(S);
    T = sizeS(2);

    F_L = dftmtx(L)/sqrt(L);

    FP = kron(F_L, eye(K));
    XT = pagemtimes(FP, S);
    X = reshape(XT, K, L, T);

elseif nargin == 3

    S = varargin{1};
    KL = varargin{2};
    K = KL(1);
    L = KL(2);
    pf = varargin{3};

    sizeS = size(S);
    T = sizeS(2);

    F_L = dftmtx(L)/sqrt(L);

    FP = kron(F_L, eye(K));
    XT = pagemtimes(FP, S);
    X = reshape(XT, K, L, T);

else
    error("Invalid number of inputs!");
end

end