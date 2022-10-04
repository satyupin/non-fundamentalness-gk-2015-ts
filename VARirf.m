function IRF = VARirf(Varest, InvA, h, Ssize)
    n = Varest.n;
    nshocks = size(InvA, 2);
    IRF = NaN(nshocks, n, h);
    for mm=1:n
        % Initialize the impulse response vector
        response = NaN(n, h);
        % Create the impulse vector
        impulse = zeros(n, 1); 
        % Set the size of the shock
        impulse(mm, 1) = Ssize; % 1 corresponds to 1 std dev shock
        
        % First period impulse response (=impulse vector)
        response(: ,1) = InvA*impulse;
        % Recursive computation of impulse response
            for kk = 2:h
                FcompN = Varest.Fcomp^(kk-1);
                response(:, kk) = FcompN(1:n, 1:n)*InvA*impulse;
            end
        IRF(mm, :, :) = response;
    end
end

