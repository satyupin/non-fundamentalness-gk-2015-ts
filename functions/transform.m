function y = transform(x, tcode)
% This function transforms the series from FRED-MD or FRED-QD databases 
% to stationary series using the allocated transformation codes.
% This function is designed specifically to be used in a loop.
%% Inputs:
% -x is an (n x 1) vector, the series to be transformed
% -tcode is a scalar with values from 1 to 7, 
%   indicating the transformation:
%   1 - No transformation
%   2 - First difference
%   3 - Second difference
%   4 - Log
%   5 - Log-difference
%   6 - Second log-difference
%   7 - First difference of percent changes (exact 6)
%% Outputs:
% -y is an (n x 1) vector, the transformed series

    N = size(x,1);
 
    % Value close to zero 
    small = 1e-6;

    % Allocate output variable
    y = NaN(N, 1);

    switch(tcode)

    case 1 % Level (i.e. no transformation): x(t)
        y = x;

    case 2 % First difference: x(t)-x(t-1)
        y(2:N) = x(2:N, 1) - x(1:N-1, 1);

    case 3 % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
        y(3:N) = x(3:N) - 2*x(2:N-1) + x(1:N-2);

    case 4 % Natural log: ln(x)
        if min(x) < small 
            y = NaN; 
        else
            y = log(x);
        end

    case 5 % First difference of natural log: ln(x)-ln(x-1)
        if min(x) > small
            x = log(x);
            y(2:N) = x(2:N) - x(1:N-1);
        end

    case 6 % Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
        if min(x) > small
            x = log(x);
            y(3:N) = x(3:N) - 2*x(2:N-1) + x(1:N-2);
        end

    case 7 % First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
        y1(2:N) = (x(2:N) - x(1:N-1))./x(1:N-1);
        y(3:N) = y1(3:N) - y1(2:N-1);

    end
end

