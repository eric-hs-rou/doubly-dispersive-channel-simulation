function [S] = OTFSmod(varargin)
% S = OTFSmod(X, pf) OTFS modulation
%   Inputs:
%      - X  : K x L x T transmit symbol tensor (K x L sized delay-Doppler grid, T symbol periods)
%      - pf : K x 1 transmit pulse-shaping filter (rectangular pulse if not provided)
%   Output:
%      - S  : KL x T matrix of OTFS modulated samples

if nargin == 1

    X = varargin{1};

    sizeX = size(X);
    K = sizeX(1);
    L = sizeX(2);
    if length(sizeX) > 2
        T = sizeX(3);
    else
        T = 1;
    end
    F_L = dftmtx(L)/sqrt(L);

    XT = reshape(X, K*L, T);
    FP = kron(F_L', eye(K));

    S = FP*XT;

elseif nargin == 2

    X = varargin{1};
    pf = varargin{2};

    sizeX = size(X);
    K = sizeX(1);
    L = sizeX(2);
    if length(sizeX) > 2
        T = sizeX(3);
    else
        T = 1;
    end

    F_L = dftmtx(L)/sqrt(L);

    XT = reshape(X, K*L, T);
    FP = kron(F_L', diag(pf));

    S = FP*XT;

else
    error("Invalid number of inputs!");
end

end