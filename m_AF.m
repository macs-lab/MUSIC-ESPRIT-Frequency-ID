function omega = m_AF(x,p)
%Anihilating Filtering: a frequency estimation method.
%-----
%USAGE	omega = m_AF(x,p)
%
%	x : input sequence
%       p : number of complex exponentials to find
%
%
%---------------------------------------------------------------
% copyright 2010, by Xu Chen. 
%---------------------------------------------------------------

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
    if omega(ii)>pi
        omega(ii) = omega(ii)-2*pi;
    end
end
