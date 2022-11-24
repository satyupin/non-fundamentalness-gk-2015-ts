function output = VARnlag(Xendo, Xexo, maxp)
% This function estimates a number of lags to be included in a VAR using 
% Akaike, Schwarz and Hannan-Quinn Information criteria
%% Inputs:
% -Xendo is an (T x n) matrix of endogeneous variables
% -Xexo is an (T x k) matrix of exogeneous variables (+ constant/trends)
% -maxp is a maximum number of lags allowed
%% Outputs:
% -ValueAIC: value of Akaike Information Criterion at the minimum
% -ValueBIC: value of Schwarz Information Criterion at the minimum
% -ValueHQ: value of Hannan-Quinn Information Criterion at the minimum
% -AICn: number of lags by Akaike Information Criterion 
% -BICn: number of lags by Schwarz Information Criterion
% -HQn: number of lags by Hannan-Quinn Information Criterion
    %% Set the data and dimensions
    
    nexo = size(Xexo, 2);
    logL = zeros(maxp, 1);
    AIC  = zeros(maxp, 1);
    BIC  = zeros(maxp, 1);
    HQ   = zeros(maxp, 1);
    for i = 1:maxp
        VARa = VARest(Xendo, Xexo, i);
        T = VARa.T;  
        n = VARa.n;
        U = VARa.U;
        % Covariance of reduced-form residuals
        Sigma = (1/(T - i*n - nexo)).*(U)'*(U);
        logL(i) = -(T*n/2)*log(2*pi)+T/2*log(det(inv(Sigma)))-(T*n/2); 
        % T*n/2 = trace(B'*inv(Sigma)*B)
        % Information criteria
        AIC(i) = log(det(Sigma))+2*(n^2*i+nexo)/T;                    
        BIC(i) = log(det(Sigma))+(n^2*i+nexo)/T*log(T);
        HQ(i)  = log(det(Sigma))+2*(n^2*i+nexo)/T*log(log(T));
    end
    [AIC, AICn] = min(AIC);
    [BIC, BICn] = min(BIC);
    [HQ, HQn] = min(HQ);
    
    %% Output
    output.ValueAIC = AIC;
    output.ValueBIC = BIC;
    output.ValueHQ  = HQ;
    output.AICn     = AICn;
    output.BICn     = BICn;
    output.HQn      = HQn;
    
end

