function output = OLSest(Y, X, NameValueArgs)
% This function estimates a regression of Y on X and returns a struct with   
% a summary of information. 
% This version is tested on R2020b since it uses namevalue args!
%% Inputs:
% -y is an (n x 1) vector, explanatory variable
% -X is an (n x k) matrix of regressors without a constant
%
% -NameValueArgs.Lag is a (1 x 1) scalar that specifies the 
%   number of NameValueArgs.Lags using in HAC standard errors computations,
%   0 default
%
% -NameValueArgs.robust is a (1 x 1) scalar that specifies the 
%   variance type in the standard errors computations: 
%   0 - standard SE (default)
%   1 - White SE
%   2 - HAC SE w/ Newey-West estimator and Bartlett kernel
% 
% -NameValueArgs.const is a (1 x 1) scalar that specifies whether the
%   constant should be included in the set of regressors
%   0 - no constant
%   1 - constant included (default)
%% Outputs:
% -beta: vector w/ OLS estimates
% -yhat: fitted values
% -resid: residuals
% -SSR: sum of squared residuals
% -SST: sum of squares total
% -SSE: sum of squares explained
% -rmse: root mean squared error
% -R2: R squared
% -R2: adj R squared
% -n: sample size
% -k: number of regressors              
% -varbeta: variance-covariance matrix of beta
% -varbetarob: robust (White or HAC) variance-covariance matrix of bhat
% -F: F statistic for test of overall significance (only
%                       well defined if there is a constant
% -Fpval: p-value of overall significance test
% -Frobust: robust F statistic (based on large sample OLS distribution)
% -Frobustpval: p-value of robust F statistic
% -lag: number of lags used in autocorrelation structure
    arguments
        Y (:, 1) {mustBeNumeric,mustBeReal}
        X (:, :) {mustBeNumeric,mustBeReal}
        NameValueArgs.const (1, 1) {mustBeNumeric,mustBeReal} = 1
        NameValueArgs.robust (1, 1) {mustBeNumeric,mustBeReal} = 0
        NameValueArgs.Lag (1, 1) {mustBeNumeric,mustBeReal} = 0
    end
    
    n = size(X, 1);
    
   

    if NameValueArgs.robust == 2 && NameValueArgs.Lag == 0
            NameValueArgs.Lag = round(4*(n/100)^(2/9));  
            % if lag number is not specified, but HAC is used,
            % set it to Newey and West rule of thumb
    end
    
    if NameValueArgs.const == 1
        X = [X ones(n,1)];
    end
    k = size(X, 2);
    beta  = inv(X'*X)*(X'*Y);
    yhat  = X*beta;
    resid = Y-yhat;
    SSR   = resid'*resid;
    SST   = (Y-mean(Y))'*(Y-mean(Y));
    SSE   = (yhat-mean(Y))'*(yhat-mean(Y));
    rmse  = sqrt(1/(n-k)*SSR);
    R2    = 1 - SSR/SST; 
    R2adj = 1 - (n-1)/(n-k)*SSR/SST;
    
    varbeta = rmse^2*inv(X'*X);
    
    switch NameValueArgs.robust
        case 1
            Shat = zeros(k);      
            for ii = 1:n
                Shat = Shat+resid(ii)^2*X(ii,:)'*X(ii,:);
            end
            varbetarob = n/(n-k)*inv(X'*X)*Shat*inv(X'*X);
        case 2
            Shat0 = zeros(k);        
            for ii = 1:n
                Shat0 = Shat0+resid(ii)^2*X(ii,:)'*X(ii,:);
            end
            Shat = Shat0;
            for l = 1:lag
                Shatl = zeros(k); 
                for ii = 1+l:n
                    Shatl = Shatl+(1-l/(lag+1))*resid(ii)*resid(ii-l)*...
                        (X(ii,:)'*X(ii-l,:)+X(ii-l,:)'*X(ii,:)); 
                    % Bartlett window
                end
                Shat = Shat + Shatl;
            end
            varbetarob = n/(n-k)*inv(X'*X)*Shat*inv(X'*X);
        otherwise
            varbetarob = [];
    end
    % F test of overall significance
    R = zeros(k-1,k);
    R(1:k-1,1:k-1) = eye(k-1);  % all coeffs zero except constant
    q = zeros(k-1,1);
    
    F = 1/(k-1)*(R*beta-q)'*inv(rmse^2*R*(inv(X'*X))*R')*(R*beta-q);         
    % small sample
    Fpval       = 1 - fcdf(F,k-1,n-k);
    
    if NameValueArgs.robust>0
        Frobust = 1/(k-1)*(R*beta-q)'*inv(n/(n-k)*R*inv(X'*X)...
            *Shat*inv(X'*X)*R')*(R*beta-q);    % large sample (Wald test)
        % Frobust = 1/(k-1)*(R*bhat-q)'*inv(R*varbetarob*R')*(R*bhat-q);    
        % large sample (Wald test)
        Frobustpval = 1 - fcdf(Frobust,k-1,n-k);
    else
        Frobust = [];
        Frobustpval = [];
    end
    output.beta    = beta;
    output.yhat    = yhat;
    output.resid   = resid;
    output.SSR     = SSR;
    output.SST     = SST;
    output.SSE     = SSE;
    output.rmse    = rmse;
    output.R2      = R2;
    output.R2adj   = R2adj;
    output.n       = n;
    output.k       = k;
    output.varbeta = varbeta;
    output.varbetarob = varbetarob;
    output.F       = F;
    output.Fpval   = Fpval;
    output.Frobust     = Frobust;
    output.Frobustpval = Frobustpval;
    output.lag     = NameValueArgs.Lag;
end

