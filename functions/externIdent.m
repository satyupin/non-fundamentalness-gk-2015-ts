function output = externIdent(U, Sigma, IV)
% This function identifies z structural shocks using external instruments
% approach. For identification the reduced-form residuals of first z
% variables from an estimated VAR are used. The estimated shock is
% normalized to induce 1 SD increase on impact.
%% Inputs:
% -U is a (n x k) matrix of reduced-form residuals from a VAR
% -Sigma is a (k x k) reduced-form covariance matrix
% -IV is a (n x z) matrix of external instruments
%% Outputs:
% -b1: a (k x z) structural impact matrix. 
% (cols of b1 are that of A0^{-1}, a matrix that post-multiplies 
% the companion form or pre-multiplies reduced-form residuals 
% in IRF computation)
% -olsEst: a struct, containing information about the first-stage
% regressions
IVn = size(IV, 2);
n = size(U, 2);
uhat = [];
for ik = 1:IVn
    inam = strcat('r',num2str(ik));
    olsEst.(inam) = OLSest(U(:,ik), IV, "const", 1, "robust", 1);

    disp('First-stage:')
    fprintf(strcat('F-stat: %4.3f, p-value: %4.3f, F-stat (robust):',...
        '%4.3f, p-value: %4.3f, R^2: %4.3f, R^2 (adj): %4.3f \n'), ...
    olsEst.(inam).F,olsEst.(inam).Fpval,olsEst.(inam).Frobust, ...
    olsEst.(inam).Frobustpval,olsEst.(inam).R2,olsEst.(inam).R2adj)
    uhat = [uhat olsEst.(inam).yhat];
end

% "Second stage" (footnote 4 in GK2015)
b21ib11_2SLS    =   [ones(length(IV),1) uhat]\U(:,IVn+1:end);  
b21ib11 = b21ib11_2SLS(2:end,:)';      % 2 SLS coefficients ignoring const               
Sig11   = Sigma(1:IVn,1:IVn);
Sig21   = Sigma(IVn+1:n,1:IVn);
Sig22   = Sigma(IVn+1:n,IVn+1:n);
ZZp     = b21ib11*Sig11*b21ib11'-(Sig21*b21ib11'+b21ib11*Sig21')+Sig22;
b12b12p = (Sig21- b21ib11*Sig11)'*(ZZp\(Sig21- b21ib11*Sig11));
b11b11p = Sig11-b12b12p;
b11     = sqrt(b11b11p);
b1      = [b11; b21ib11*b11];
% b1unit  = [1; b21ib11]*shockSize;

output.b1 = b1;
output.Fstage = olsEst;
end

