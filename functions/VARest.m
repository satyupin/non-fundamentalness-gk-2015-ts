function output = VARest(Xendo, Xexo, p)
% This function estimates a VAR(p) model and returns a struct with a  
% summary of information
%% Inputs:
% -Xendo is an (T x n) matrix of endogeneous variables
% -Xexo is an (T x k) matrix of exogeneous variables (+ constant/trends)
% -p is a number of lags specified by the user
%% Outputs:
% -Y: matrix of all left-hand side variables (in the order specified in
% Xendo)
% -X: matrix of all right-hand side variables, exo and all lags of endo
% -Xexo: matrix of all exogeneous variables (order preserved)
% -B: matrix of all OLS coefficients (in the order of X)
% -Fcomp: companion form matrix (used to compute IRFs and eigenvalues)
% -Yhat: fitted values
% -maxEig: the highest eigenvalue (stability check)
% -U: residuals (in the order of Y)
% -Sigma: reduced-form covariance matrix
% -VarB: covariance matrix of coefficients (heteroskedasticity assumed)
% -logL: log-likelihood
% -AIC: value of Akaike Information Criterion
% -BIC: value of Schwarz Information Criterion
% -HQ: value of Hannan-Quinn Information Criterion
%                       well defined if there is a constant
% -T: original sample size
% -n: number of endogenous variables
% -p: number of lags (specified by the user)
% -nexo: number of exogenous variables per equation
    %% Set the data and dimensions
    [T, n] = size(Xendo);
    % Create lags
    X = lagmatrix(Xendo,1:p); % define right hand variables
    Xorig = X;
    X = X(p+1:end,:); 
    XexoO = Xexo;
    Xexo = Xexo(p+1:end,:);
    nexo = size(Xexo, 2);
    Yorig = Xendo;
    X = [Xexo X];
    Y = Xendo(p+1:end, :);

    %% Estimate coefficients and IC

    B    = (X'*X)\(X'*Y); % OLS 
    Yhat = X*B;         % fitted values
    % Companion form
    Bt = B';
    Fcomp = [Bt(:,1+nexo:n*p+nexo); eye(n*(p-1)) zeros(n*(p-1),n)]; 
    maxEig = max(abs(eig(Fcomp))); % stability
    U = Y-Yhat; % residuals
    Sigma = U'*U/(T-p*n-nexo); % covariance matrix
    VarB = kron(inv(X'*X),Sigma); % covariance matrix of coefs
    % Log-likelihood
    logL = -(T*n/2)*log(2*pi)+T/2*log(det(inv(Sigma)))-(T*n/2); 
    % T*n/2 = trace(B'*inv(Sigma)*B)
    % Information criteria
    AIC = log(det(Sigma))+2*(n^2*p+nexo)/T;                    
    BIC = log(det(Sigma))+(n^2*p+nexo)/T*log(T);
    HQ  = log(det(Sigma))+2*(n^2*p+nexo)/T*log(log(T));
    
    %% Output
    output.Y     = Y;
    output.X     = X;
    output.Xexo  = Xexo;
    output.Yorig = Yorig;
    output.Xorig = Xorig;
    output.XexoO = XexoO;
    output.B     = B';
    output.Fcomp = Fcomp;
    output.Yhat  = Yhat;
    output.maxEig = maxEig;
    output.U     = U;
    output.Sigma = Sigma;
    output.VarB  = VarB;
    output.logL  = logL;
    output.AIC   = AIC;
    output.BIC   = BIC;
    output.HQ    = HQ;
    output.T     = T;
    output.n     = n;
    output.p     = p;
    output.nexo  = nexo;
end

