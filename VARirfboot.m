function output = VARirfboot(EstVAR, IV, IRF, nsim, h, alpha, IVStartInd,...
    IVEndInd, SampStartInd)

T_ef = EstVAR.T - EstVAR.p;
Tiv = size(IV, 1);
n = EstVAR.n;
nexo = EstVAR.nexo;
IVCount = zeros(nsim, size(IV, 2));
bootIRFs = NaN(n, n, h, nsim);
p = EstVAR.p;
for j = 1:nsim
    
    index = randsample(1:T_ef, T_ef, true)';    
    bootU = EstVAR.U(index, :);
        
    bootY = zeros(EstVAR.T, EstVAR.n); 
    bootY(1:p, :) = EstVAR.Yorig(1:p, :); %initial values of y, same for all j
    for i = p+1:T_ef+p
        bootY(i, :)= EstVAR.B*[EstVAR.Xexo(i-p, :)'; ...
            vec(fliplr(bootY(i-p:i-1, :)'))] + bootU(i-p, :)'; % bootstrap
    end
    bootvarEst = VARest(bootY, EstVAR.XexoO, EstVAR.p);
        
    % get residuals for identification sample
    bootUhat = bootvarEst.U(IVStartInd-SampStartInd-p:...
        IVEndInd-SampStartInd-p,:);
    
    index2 = randsample(1:Tiv,Tiv,true)';
    bootIV = IV(index2, :);
    bootU = bootUhat(index2, :);
    bootSigma = bootU'*bootU/(EstVAR.T-p*n-nexo);

    % count the number proxy variables not censored to zero
    IVCount(j, :) = sum(abs(bootIV) > 0,1);

    if IVCount(j, :)<15
        continue
    end
    bootExIV = externIdent(bootU, bootSigma, bootIV);
    bootb1 = bootExIV.b1;
    bootInvA = zeros(n, n);
    bootInvA(:, 1) = bootb1;
    % compute IRFs
    bootIRFs(:, :, :, j)  = VARirf(bootvarEst, bootInvA, h, 1);
end
IRFsmed = quantile(bootIRFs, 0.5, 4);

IRFsupper = quantile(bootIRFs, 1-alpha/2, 4);  
IRFslower = quantile(bootIRFs, alpha/2, 4);

output.IRFmed = IRFsmed;
output.IRFmean = IRF;
output.IRFupper = IRFsupper;
output.IRFlower = IRFslower;
end

