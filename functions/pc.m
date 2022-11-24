function [ehat, fhat, lambda, ss] = pc(Y, nfac)
% principal components with normalization F'F/T=I
% Y is observed (n x k) matrix containing series
% r is the true number of true factors
% F is T by r matrix of true factors
% Lambda N by r is the true loading matrix
% C=F*Lambda' T by N is the true common component
% chat is the estimated common component
% The function is taken from Fabio Ferroni's Github page (BVAR toolbox)
[bigt, bign] = size(Y);
[Fhat0, eigval, ~] = svd(Y*Y');
fhat = Fhat0(:, 1:nfac)*sqrt(bigt);
lambda = Y'*fhat/bigt;

%chi2=fhat*lambda';
%diag(lambda'*lambda)
%diag(fhat'*fhat)                % this should equal the largest eigenvalues
%sum(diag(eigval(1:nfac,1:nfac)))/sum(diag(eigval))
%mean(var(chi2))                 % this should equal decomposition of variance

ehat = Y - fhat*lambda';

ve2 = sum(ehat'.*ehat')'/bign;
ss = diag(eigval);
end

