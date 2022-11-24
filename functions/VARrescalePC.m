function [ReScale] = VARrescalePC(Lambda, n1, orderPC, Scale)
% This function obtains a rescaling factor for principal components
% used in a FAVAR estimation to revert the IRFs back to original variables.
% The function is taken from Fabio Ferroni's Github page (BVAR toolbox).
%% Inputs:
% -Lambda is a (n2 x nFac) matrix of factor loadings
% -n1 is a scalar, a number of variables included in FAVAR not as factors
% -orderPC is a scalar, taking values of:
%   1 - factors are ordered first in FAVAR estimation
%   2 - factors are ordered second in FAVAR estimation
% -Scale is an optional (n2 x 1) vector of standard deviations of the 
% original vars in case they weren't rescaled before calling this function
%% Outputs:
% -Rescale is a scaling factor to be pre-multiplied with the IRFs to obtain
% responses of the original variables
if nargin < 4
    Scale = ones(size(Lambda, 1), 1);
end
nFac  = size(Lambda,2);
n2  = size(Lambda,1);

ReScale = NaN;
    
    
switch orderPC
    case 1 % factor first

        Scale_  = repmat([Scale; ones(n1,1)], 1, nFac + n1);
        Lambda_ = [Lambda zeros(n2 , n1); zeros(n1, nFac + n1)];
        Lambda_(n2+1:n2+n1, nFac+1:nFac+n1) = eye(n1);
        ReScale = Scale_ .* Lambda_ ;

    case 2 % factor second
        
        Scale_  = repmat([ones(n1, 1); Scale], 1, nFac + n1);
        Lambda_ = [zeros(n1, n1+nFac); zeros(n2 , n1) Lambda];
        Lambda_(1:n1, 1:n1) = eye(n1);
        ReScale = Scale_ .* Lambda_ ;

end
