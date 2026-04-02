function [omega, sigma] = AF_with_denoising(x,p,L,maxIter)
%Anihiliting Filter with Catzow's denoising for Frequency estimation.
%-----
%USAGE	Px=AF_with_denoising(x,p,L)
%
%	The input sequence x is assumed to consist of p complex
%	exponentials in white noise.  The frequencies of the
%	complex exponentials and the variance of the white noise
%	are estimated using the MUSIC algorithm.
%
%	x : input sequence
%       p : number of complex exponentials to find
%	L : number of columns for Toeplitz matrix
%
%
%  see also PHD, EV, and MIN_NORM
%
%---------------------------------------------------------------
% copyright 2010, by Xu Chen. 
%---------------------------------------------------------------

x   = x(:);
% x = flipud(x);
N   = length(x);
if nargin == 3
    maxIter = 100;
end
if L<p+2 || N<L, error('Size of R is inappropriate'), end
x   = x - ones(N,1)*(sum(x)/N);
Y   = toeplitz( x(L:end) , flipud( x(1:L) ) ); % Toeplitz matrix
[U,S,V]     = svd(Y); % Y = U*S*V'
sValue      = diag(S);
sValue      = sValue(sValue~=0); % get the singular values
sigma       = zeros(L,maxIter); % expected to perform less than 100 iterations in denoising
sigma(:,1)  = sValue;
% figure;semilogy(sValue,'*-');
% xlabel('Index on the diagonal');ylabel('Singular values \sigma')
% title('Original data');

sDiff       = diff(sValue);
[~,sDiffMaxPos] = max(-sDiff); % abs(sDiff)-->-sDiff
if sDiffMaxPos<p
    disp('Number of narrow band disturbances less than specified, reduction performed in Catzow denoising');
    p = sDiffMaxPos
end
iter        = 0;
while(1)
    if -sDiff(p)>-80*sDiff(p+1) || iter >= maxIter % recall that sigma1>sigma2>...
%         abs(sDiff(p))>5*abs(sDiff(p+1))
%         sValue(p)/sValue(p+1)>2*sValue(p-1)/sValue(p) % may need a better threshold
        break;
        % S_hat = zeros(size(S));
        % for ii = 1:p
        %     S_hat(ii,ii) = S(ii,ii);
        % end
    else
        iter    = iter + 1;
        sizeS   = size(S);
        Y_hat   = zeros(sizeS);
        % Y_hat = zeros(N-L,L);
        for ii = 1:p
            Y_hat = Y_hat + S(ii,ii)*U(:,ii)*V(:,ii)'; % principle component
        end
        % Averaging to get the Toeplitz form
        for ii = sizeS(2)-1:-1:-sizeS(1)+1
            x(ii+sizeS(1)) = mean(diag(Y_hat,ii));
        end
        Y_hat_toeplitz = toeplitz(...
            x(L:end),flipud( x(1:L) )...
            ); % Toeplitz matrix
        [U,S,V] = svd(Y_hat_toeplitz); % Y = U*S*V'
        sValue  = diag(S);
        sValue  = sValue(sValue~=0); % get the singular values
        sDiff   = diff(sValue);
        sigma(:,iter+1) = sValue;
    end
end

sigma(:,iter+2:end) = []; % erase the unused columns
sprintf('Number of iterations performed in Catzow denoising: %d',iter)
%% anihilating filtering
Y_AF        = toeplitz(...
    x(p+1:end),flipud(x(1:p+1))); % Toeplitz matrix
[~,~,V_AF]  = svd(Y_AF); % Y = U*S*V'
h           = V_AF(:,p+1);
zos         = roots(h);
omega       = -phase(zos);
for ii = 1:p
%     if abs(omega(ii))>pi
%         omega(ii) = abs(omega(ii))-pi;
%     end
    while(omega(ii)<0)
        omega(ii) = omega(ii)+2*pi;
    end
end
