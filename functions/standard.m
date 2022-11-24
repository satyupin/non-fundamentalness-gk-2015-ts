function [X, Scale, Loc] = standard(Y)
% This function transforms the series by subtracting the mean and dividing
% by the standard deviation.
%% Inputs:
% -Y is an (n x k) matrix, the columns of which contain the series 
% to be standardized
%% Outputs:
% -X is an (n x k) matrix, the columns of which contain the transformed
% series
% -Scale is a (k x 1) vector of standard deviations of the original series
% -Loc is a (k x 1) vector of means of the original series
    T = size(Y ,1);
    my = repmat(nanmean(Y), T, 1);
    sy = repmat(nanstd(Y), T, 1);
    X = (Y - my)./sy;
    Scale = nanstd(Y)';
    Loc = nanmean(Y)';
end

