function [x_denoised, sigma] = Denoising(x,p,L,maxIter,gain,selfjudge)
%Cadzow's denoising method. (Essentially a principle component method)
%-----
%USAGE	[x_denoised, sigma] = Denoising(x,p,L,maxIter)
%
%	The input sequence x is assumed to consist of p complex
%	exponentials in white noise.  The frequencies of the
%	complex exponentials and the variance of the white noise
%	are estimated using the MUSIC algorithm.
%
%	x : input sequence
%       p : number of complex exponentials to find
%	L : number of columns for Toeplitz matrix
%   maxIter: maximum iteration allowed
%   gain: Gain of seperation of noise sigular values from signal singular
%   values
%
%
%  see UC Berkeley EE225A Notes
%  see also PHD, EV, and MIN_NORM
%
%---------------------------------------------------------------
% copyright 2010, by Xu Chen. 
%---------------------------------------------------------------

if nargin   <= 3
    maxIter = 100;
    gain = 80;
elseif nargin<=4
    gain = 80;
end
x           = x(:);
% x = flipud(x);
N           = length(x);
if L<p+2 || N<L, error('Size of R is inappropriate'), end
% x           = x - ones(N,1)*(sum(x)/N);
Y           = toeplitz( x(L:end) , flipud( x(1:L) ) ); % Toeplitz matrix
[U,S,V]     = svd(Y); % Y = U*S*V'
sValue      = diag(S);
sValue      = sValue(sValue~=0); % get the singular values
sigma       = zeros(L,maxIter); % expected to perform less than 100 iterations in denoising
sigma(:,1)  = sValue;
% figure;semilogy(sValue,'*-');
% xlabel('Index on the diagonal');ylabel('Singular values \sigma')
% title('Original data');

sDiff           = diff(sValue);
[~,sDiffMaxPos] = max(-sDiff); % abs(sDiff)-->-sDiff
if strcmp(selfjudge,'selfjudge')
    if sDiffMaxPos<p
        disp('Number of narrow band disturbances less than specified, reduction performed in Catzow denoising');
        p           = sDiffMaxPos
    end
end
iter = 0;
while(1)
    if -sDiff(p)>-gain*sDiff(p+1) || iter >= maxIter % recall that sigma1>sigma2>...
        break;
    else
        iter        = iter + 1;
        sizeS       = size(S);
        Y_hat       = zeros(sizeS);
        for ii      = 1:p
            Y_hat   = Y_hat + S(ii,ii)*U(:,ii)*V(:,ii)'; % principle component
        end
        % Averaging to get the Toeplitz form
        for ii      = sizeS(2)-1:-1:-sizeS(1)+1
            x(ii+sizeS(1)) = mean(diag(Y_hat,ii));
        end
        Y_hat_toeplitz  = toeplitz(...
            x(L:end),flipud( x(1:L) )...
            ); % Toeplitz matrix
        [U,S,V]         = svd(Y_hat_toeplitz); % Y = U*S*V'
        sValue          = diag(S);
        sValue          = sValue(sValue~=0); % get the singular values
        sDiff           = diff(sValue);
        sigma(:,iter+1) = sValue;
    end
end
x_denoised              = x;
sigma(:,iter+2:end)     = []; % erase the unused columns

sprintf('Number of iterations performed in Catzow denoising: %d',iter)
