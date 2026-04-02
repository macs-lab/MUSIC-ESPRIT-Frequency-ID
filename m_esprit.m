function omega = m_esprit(x,p,M)
%ESPRIT	Frequency estimation using the ESPRIT algorithm.
%-----
%USAGE	Px=m_esprit(x,p,M)
%
%	The input sequence x is assumed to consist of p complex
%	exponentials in white noise.  The frequencies of the
%	complex exponentials and the variance of the white noise
%	are estimated using the MUSIC algorithm.
%
%	x : input sequence
%       p : number of complex exponentials to find
%	M : number of noise eigenvectors to use
%
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
% v(:,i(j) (for j = M-p+1:M) is the signal eigenvector

for ii = 1:p
    if y(M-ii)/y(M-ii+1)<0.2
        p = ii;
        disp('Number of narrow-band signals less than its assigned value. Reduction was performed.')
        break
    end
end
Qs          = zeros(M,p);
Qs_upbar    = zeros(M-1,p);
Qs_downbar  = zeros(M-1,p);
for ii = 1:p
    Qs(:,ii)            = v(:,i(M-ii+1));
    Qs_upbar(:,ii)      = v(2:end,i(M-ii+1));
    Qs_downbar(:,ii)    = v(1:end-1,i(M-ii+1));
end
phi         = (Qs_downbar'*Qs_downbar)\Qs_downbar'*Qs_upbar;
% phi = inv(Qs_downbar'*Qs_downbar)*Qs_downbar*Qs_upbar;
lambda      = eig(phi);
% omega = -log(lambda)./1i; % use this if data flipping was applied at line 25
omega       = log(lambda)./1i; % use this if no data flipping was applied at line 25
omega       = real(omega);
for ii = 1:p
    if omega(ii)<0
        omega(ii) = omega(ii)+2*pi;
    end
end
