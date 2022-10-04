function [Data, vnames, tcodes, addinfo] = removeVars(DataStruct,...
    VarsToRemove, vnames, tcodes, addinfo)
    if nargin < 5
        VarInfo = [vnames; num2cell(tcodes)];
    end
    VarInfo = [vnames; addinfo; num2cell(tcodes)];
    
    
    [~,indx] = ismember(VarsToRemove,VarInfo(1, :));
    VarInfo(:, indx) = [];
    
    vnames = VarInfo(1, :);
    addinfo = VarInfo(2:end-1, :);
    tcodes = cell2mat(VarInfo(end, :));
    for i = 1:length(vnames)
        Data.(vnames{i}) = DataStruct.(vnames{i});
    end
end

