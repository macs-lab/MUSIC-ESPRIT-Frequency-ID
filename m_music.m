function Px = m_music(x,p,M,omegaVector)
%MUSIC	Frequency estimation using the MUSIC algorithm.
%-----
%USAGE	Px=music(x,p,M)
%
%	The input sequence x is assumed to consist of p complex
%	exponentials in white noise.  The frequencies of the
%	complex exponentials and the variance of the white noise
%	are estimated using the MUSIC algorithm.
%
%	x : input sequence
%       p : number of complex exponentials to find
%	M : number of noise eigenvectors to use
%   omegaVector: a priori information of the narrow band frequency region
%
%	The frequency estimates are found from the peaks of the
%	pseudospectrum Px.
%
%  see also PHD, EV, and MIN_NORM
%
%---------------------------------------------------------------
% copyright 2010, by Xu Chen. 
%---------------------------------------------------------------

x           = x(:);
% x = flipud(x);
N           = length(x);
if M<p+1 || N<M, error('Size of R is inappropriate'), end
x           = x - ones(N,1)*(sum(x)/N);
R           = convm(x,M)'*convm(x,M)/(N-1);
[v,d]       = eig(R); % v is the eigenvector matrix; d is the diagonal eigenvalue matrix
[y,i]       = sort(diag(d));
% e.g.
% d =
%    1.000000000000000                   0                   0
%                    0   2.000000000000000                   0
%                    0                   0   0.100000000000000
% [y,i] = sort(diag(d))
% y =
%    0.100000000000000
%    1.000000000000000
%    2.000000000000000
% i =
%      3
%      1
%      2
%
% v(:,i(j) (for j = 1:M-p) is the noise eigenvector
Px          = zeros(length(omegaVector),1);
a           = zeros(M,1);
for jj = 1:length(omegaVector)
    for ii = 1:M
        a(ii) = exp(1i*(ii-1)*omegaVector(jj)); % use this if x is not flipped in line 26
%         a(ii) = exp(-1i*(ii-1)*omegaVector(jj)); % use this if x is
%         flipped in line 26.
    end
    for kk = 1:M-p
        Px(jj) = Px(jj) + abs(v(:,i(kk))'*a)^2;
    end
end
Px=-20*log10(Px);
% Px=1./Px;
