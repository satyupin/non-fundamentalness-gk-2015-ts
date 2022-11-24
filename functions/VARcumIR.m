function IRF = VARcumIR(IRF, cumid)
    for j = 1:size(IRF, 1)
        for i = 1:size(IRF, 2) % compute CIRF if requested
            if cumid(i) == 1
                IRF(j, i, :) = IRF(j, i, :);
            elseif cumid(i) == 2
                IRF(j, i, :) = cumsum(IRF(j, i, :));
            elseif cumid(i) == 3
                IRF(j, i, :) = cumsum(cumsum(IRF(j, i, :)));
            elseif cumid(i) == 4
                IRF(j, i, :) = exp(IRF(j, i, :));
            elseif cumid(i) == 5
                IRF(j, i, :) = exp(cumsum(IRF(j, i, :)));
            elseif cumid(i) == 6 || cum(i) == 7
                IRF(j, i, :) = exp(cumsum(cumsum(IRF(j, i, :))));
            end
        end
    end
end

