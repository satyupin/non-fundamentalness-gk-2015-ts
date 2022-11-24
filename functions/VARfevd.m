function FEVD = VARfevd(Varest, h, InvA)

n = Varest.n;
p = Varest.p;
nshocks = n;
Sigma = Varest.Sigma; %or do I use the truncated sample?
F = Varest.Fcomp;
VARk = [Sigma zeros(n,n*(p-1)); zeros(n*(p-1),n) zeros(n*(p-1),n*(p-1))];


MSE(:, :, 1) = VARk;

for kk = 2:h
   VARk = F*VARk*F';
   MSE(:, :, kk) = MSE(:, :, kk-1) + VARk;
end

FEVD = zeros(n, nshocks, h);

for j = 1:nshocks

VARk_j = [InvA(:, j)*InvA(:, j)' zeros(n, n*(p-1)); ...
    zeros(n*(p-1), n) zeros(n*(p-1), n*(p-1))];
MSE_j(:, :, 1) = VARk_j;

for kk = 2:h
    VARk_j = F*VARk_j*F';
    MSE_j(:, :, kk) = MSE_j(:, : ,kk-1) + VARk_j;   
end
% Forecast Error Covariance Decomp
FECD = MSE_j./MSE;

for nn = 1:h
   FEVD(:, j, nn) = diag(FECD(1:n, 1:n, nn));  
end
end
end

