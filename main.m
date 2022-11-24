
clear all
close all
clc

addpath(genpath('functions'))
addpath(genpath('data'))
%% Download original GK2015 data
[xlsdata, xlstext] = xlsread('GK2015_Data.xlsx','Sheet1');
dates = xlstext(3:end, 1);
datesnum = Date2Num(dates, 'm');
% Variable names
vnames_long = xlstext(1, 2:end);
vnames = xlstext(2, 2:end);
nvar = length(vnames);
% Remove 123456789 to NaN
data   = Num2NaN(xlsdata);
for ii=1:length(vnames)
    DATA.(vnames{ii}) = data(:, ii);
end

%% Download FRED-MD
md = importdata('fred_md_012022.csv', ',');
mdNames = importdata('FRED-MD_updated_appendix.csv', ',');
% Variable names
vnamesMdLong = (mdNames.textdata(2:end, 4))';
vnamesMd = md.textdata(1, 2:end);
% Transformation codes
tcodes = md.data(1, :);

% Download the data
mdData = md.data(2:end, :);

% Convert dates to the same format
mdDates = cell(size(mdData, 1), 1);
for dd = 1:size(mdData, 1)
    tempDatevec = datevec(md.textdata(dd+2, 1));
    tempDate = ConvDate(tempDatevec);
    mdDates{dd} = tempDate;
end
% Variables info in one array 
VarInfo = [vnamesMd; vnamesMdLong; num2cell(tcodes)];

%% Set the parameters
p = 12; % lags
h = 40; % IRF horizon
alpha = 0.1; % confidence bands
smplStart = dates{1};
mdStart = '1978m2';
smplEnd = dates{end};
mdEnd = '2021m6';
IVsmplStart = '1991m1'; % IV length (full sample is used for estimation)
IVsmplEnd = smplEnd;
mdStartInd = find(strcmp(mdDates, mdStart));
mdEndInd = find(strcmp(mdDates, mdEnd));
TT = mdEndInd - mdStartInd + 1; % Effecitve size of FRED-MD
NN = size(mdData, 2);
T = size(data,1); % Size of the original GK2015 dataset (for replication)
nFac = 11; % factors
nrep = 500; 
computeBT = true; % bootstrap
removeTrend = true;
dataFrequency = 'M';
%% Transform and save the data
% Populate the struct with stationary series
for ii = 1:length(vnamesMd)
    temp = transform(mdData(:, ii), tcodes(ii));
    temp = fillmissing(temp, 'previous');
    temp = temp(mdStartInd:mdEndInd);
    
    DATAMD.(vnamesMd{ii}) = temp;
end
%% Set variable matrices

% Set endogenous
VARvnames_long = {'1yr rate';'CPI';'IP';'EBP';};
VARvnames      = {'GS1';'CPIAUCSL';'INDPRO';'EBP';};
n       = length(VARvnames);

% Set IV
IVvnames_long = {'FF4';};
IVvnames      = {'FF4'};
IVn        = length(IVvnames);

% Set exogeneous
const = ones(TT, 1);
Xexo = [const];
nexo = size(Xexo, 2);

% Create matrices of variables for the VAR
Xendo = nan(TT,n);
for ii = 1:n
    Xendo(:, ii) = DATAMD.(VARvnames{ii});
end

% Create matrices of variables for the instrument
IV = nan(TT,IVn);
for ii = 1:IVn
    IV(:, ii) = DATAMD.(IVvnames{ii});
end

% Align the sample
%smplStartInd = find(strcmp(dates, smplStart));
%smplEndInd = find(strcmp(dates, smplEnd));

IVStartInd = find(strcmp(mdDates, IVsmplStart));
IVEndInd = find(strcmp(mdDates, IVsmplEnd));
IV = IV(IVStartInd-mdStartInd:IVEndInd-mdStartInd, :);

% CPI and IP are in logs in GK2015
%Xendo(:, 2:3) = log(Xendo(:, 2:3));
%% Make PC
% Prepare data 
% here I remove series that are either: 
%   -present in a VAR as they are (GS1, INDPRO, EBP, CPI)
%   -contain info for identification and are assumed to be external (FF4)
%   -have a lot of missing data (key results do not change much if include)
[DATAPC, vnamesPC, tcodesPC, vnamesPCLong] = removeVars(DATAMD, ...
    {'GS1'; 'FF4'; 'NONBORRES'; 'ACOGNO'; 'CPIAUCSL'; 'INDPRO';'EBP';},...
    vnamesMd, tcodes, vnamesMdLong);
if removeTrend
    for ii = 1:length(vnamesPC)
            temp = DATAPC.(vnamesPC{ii});
            tr = [1:size(temp, 1)]';
            trEst = OLSest(temp, tr);
            trEsts.(vnamesPC{ii}) = trEst;
            temp = trEst.resid;
            DATAPC.(vnamesPC{ii}) = temp;
    end
end
% Location and scale from standardization
PCSeries = NaN(TT, length(vnamesPC));
Scale = NaN(length(vnamesPC), 1);
Loc = NaN(length(vnamesPC), 1);
% Get standardized series
for jj = 1:length(vnamesPC)
    [PCSeries(:, jj), Scale(jj), Loc(jj)]  = ...
        standard(DATAPC.(vnamesPC{jj}));
end
% Make principal components
[~, fhat, Lambda, ~] = pc(PCSeries, nFac);

%% Estimate VAR
Varest = VARest(Xendo, Xexo, p);

% Align the sample to conform with IV
U = Varest.U(IVStartInd-p-mdStartInd:IVEndInd-p-mdStartInd,:);      
% Estimate new Sigma
Sigma = U'*U/(T-p*n-nexo);

% use proxy sample for structural shock identification 
% (potentially a subset of the estimation sample)

%% Use proxy for identification
% Compute b1 as in GK2015, footnote 4
ExIV = externIdent(U, Sigma, IV);
b1 = ExIV.b1;
InvA    = zeros(n, n);
InvA(:, 1) = b1;

% Compute struct. shock series (as in Stock and Watson (2018), ftnote 6)
StrU = (b1'*inv(Sigma)*U')'*1;

% compute IRFs to shock

IRFs_gk = VARirf(Varest, InvA, h, 1);
% cumulate them
[~,cumid] = ismember(VARvnames,vnamesMd);
cum = tcodes(cumid);
IRFs_gk = VARcumIR(IRFs_gk, cum);
%IRFs_gk(:, 2:3, :) = IRFs_gk(:, 2:3, :)/100;
IRFs_single = squeeze(IRFs_gk(1, :, :))';
if computeBT
IRF_gk_bt = VARirfboot(Varest, IV, IRFs_gk, nrep, h, alpha, IVStartInd,...
    IVEndInd, mdStartInd);
end
IRFs_singleUp = squeeze(IRF_gk_bt.IRFupper(1, :, :))';
IRFs_singleUp(:, 2:3) = (IRFs_singleUp(:, 2:3))*100;
IRFs_singleLo = squeeze(IRF_gk_bt.IRFlower(1, :, :))';
IRFs_singleLo(:, 2:3) = (IRFs_singleLo(:, 2:3))*100;
IRFs_singleMed = squeeze(IRF_gk_bt.IRFmed(1, :, :))';
IRFs_singleMed(:, 2:3) = (IRFs_singleMed(:, 2:3))*100;

FEVD_gk = VARfevd(Varest, h, InvA);
% IP and CPI are in logs in GK2015, so that measured IRFs are in %
IRFs_single(:, 2:3) = log(IRFs_single(:, 2:3))*100;


%% Orthogonality tests
% Initialize the matrix w/ results
OrthPvalues = NaN(4, nFac-1);
ii = 0; % the counter
for i = [1, 4, 8, 12] % number of lags of PCs to include in the regression
    ii = ii+1;
    for j = 1:nFac-1
        FhatL = lagmatrix(fhat(:, 1:j), 1:i);
        FhatL = FhatL(IVStartInd-p-mdStartInd:IVEndInd-p-mdStartInd, :);
        Fp = OLSest(StrU, FhatL, "const", 1, "robust", 1);
        OrthPvalues(ii, j) = Fp.Frobustpval;
    end
end


%% Number of factors in a FAVAR

OrthPvaluesF = NaN(nFac-1, nFac-1);

for i = 0:nFac-2
    if i == 0
        XendoF = Xendo;
    else
        XendoF = [Xendo(:, 1:4) fhat(:, 1:i)];
    end
    % nlagIC = VARnlag(XendoF, Xexo, p);
    NewP = 12;
    Favar = VARest(XendoF, Xexo, NewP);
    % Align the sample to conform with IV
    NewU = Favar.U(IVStartInd-NewP-mdStartInd:IVEndInd-NewP-mdStartInd,:);      
    % Estimate new Sigma
    SigmaF = NewU'*NewU/(T-NewP*(n+i)-nexo);
    % New structural shock series
    ExIVf = externIdent(NewU, SigmaF, IV);
    b1F = ExIVf.b1;
    NewStrU = (b1F'*inv(SigmaF)*NewU')'*1;
    for j = i+1:nFac-1
        % Add PCs
        FhatL = lagmatrix(fhat(:, i+1:j), 1);
        FhatL = ...
            FhatL(IVStartInd-NewP-mdStartInd:IVEndInd-NewP-mdStartInd, :);
        Favp = OLSest(NewStrU, FhatL, "const", 1, "robust", 1);
        % Save p-value from a robust F-test
        OrthPvaluesF(i+1, j) = Favp.Frobustpval;
    end
end

%% FAVAR estimation

FAVARpc = 4;
XendoF = [Xendo(:, 1:4) fhat(:, 1:FAVARpc)];
% Rescale the series back to original
OrigSeries = PCSeries.*repmat(Scale', size(PCSeries, 1), 1)...
    + repmat(Loc', size(PCSeries, 1), 1);
if removeTrend
    for ii = 1:length(vnamesPC)
    OrigSeries(:, ii) = OrigSeries(:, ii)+[ones(TT, 1) [1:TT]']*...
        trEsts.(vnamesPC{ii}).beta;
    end
end
% Obtain projection coefficients
Lhat = (fhat(:, 1:FAVARpc)'*fhat(:, 1:FAVARpc))\...
    (fhat(:, 1:FAVARpc)'*OrigSeries);
% nlagIC = VARnlag(XendoF, Xexo, p);
NewP = 12;
Favar = VARest(XendoF, Xexo, NewP);
% Align the sample to conform with IV
NewU = Favar.U(IVStartInd-NewP-mdStartInd:IVEndInd-NewP-mdStartInd,:);      
% Estimate new Sigma
SigmaF = NewU'*NewU/(T-NewP*(n+nFac)-nexo);
ExIVf = externIdent(NewU, SigmaF, IV);
b1F = ExIVf.b1;
NewStrU = (b1F'*inv(SigmaF)*NewU')'*1;
InvAf    = zeros(Favar.n, Favar.n);
InvAf(:, 1) = b1F;


%% Compute FAVAR IRFs

% Compute plain IRFs
IRFsF = VARirf(Favar, InvAf, h, 1);
% Obtain scaling factor to get the original variables from PCs
ReScale = VARrescalePC(Lhat', Favar.n-FAVARpc, 2);
IRFsOr = NaN(size(ReScale, 2), size(ReScale, 1), h);
    for k = 1:size(IRFsF, 1) % number of shocks
        if computeBT
            IRFsBTf = VARirfboot(Favar, IV, IRFsF, nrep, h, alpha, ...
            IVStartInd, IVEndInd, mdStartInd);
            IRFsOr(k, :, :) = ReScale*squeeze(IRFsBTf.IRFmean(k, :, :));
        else
            IRFsOr(k, :, :) = ReScale*squeeze(IRFsF(k, :, :));
        end
    end
% Specify variables of interest
vInt = {'COMPAPFFx'; 'GS5'; 'GS10'};
% Obtain their ids
[~,varid] = ismember(vInt, vnamesPC);
cumid = tcodesPC(varid); % [1, 5, 4, 1]
cumid = [[1, 4, 4, 1] cumid];
varid = [1; 2; 3; 4; varid];
% Cumulate IRFs to get the responses to original variables
IRFsOr = VARcumIR(IRFsOr(:, varid, :), cumid);
% Get the responses for MP shock (we only have them anyway)
IRFSingle = squeeze(IRFsOr(1, :, :))';

% Make all the variables in % (take logs of CPI and IP)
IRFSingle(:, 2:3) = log(IRFSingle(:, 2:3))*100;

FEVD_F = VARfevd(Favar, h, InvAf);

%% FEVD comparison (w/o bootstrapping)
% compare variance of output explained
ind = [10, 20, 30, 40];

FEVD_Fc = squeeze(FEVD_F(3, 1, :));
FEVD_gkc = squeeze(FEVD_gk(3, 1, :));

FEVD_comp = [FEVD_Fc  FEVD_gkc];
FEVD_comp = FEVD_comp(ind, :);
%% Compare specifications
% IRFSingle3 = IRFSingle;
IRFSingle4 = IRFSingle;


IRFSingle4Up = squeeze(IRFsBTf.IRFupper(1, :, :))';
IRFSingle4Up = (ReScale*IRFSingle4Up')';
IRFSingle4Up = VARcumIR(IRFSingle4Up(:, varid, :), cumid);
IRFSingle4Up(:, 2:3) = log(IRFSingle4Up(:, 2:3))*100;

IRFSingle4Lo = squeeze(IRFsBTf.IRFlower(1, :, :))';
IRFSingle4Lo = (ReScale*IRFSingle4Lo')';
IRFSingle4Lo = VARcumIR(IRFSingle4Lo(:, varid, :), cumid);
IRFSingle4Lo(:, 2:3) = log(IRFSingle4Lo(:, 2:3))*100;

% Bias correction (center around estimates)

IRFSingle4Med = squeeze(IRFsBTf.IRFmed(1, :, :))';
IRFSingle4Med = (ReScale*IRFSingle4Med')';
IRFSingle4Med = VARcumIR(IRFSingle4Med(:, varid, :), cumid);
IRFSingle4Med(:, 2:3) = log(IRFSingle4Med(:, 2:3))*100;

IRFSingle4Lo = IRFSingle4Lo-IRFSingle4Med+IRFSingle4;
IRFSingle4Up = IRFSingle4Up-IRFSingle4Med+IRFSingle4;

IRFs_singleLo = IRFs_singleLo-IRFs_singleMed+IRFs_single;
IRFs_singleUp = IRFs_singleUp-IRFs_singleMed+IRFs_single;

TitlenamesF = {'1yr rate';'CPI';'IP';'EBP'; 'Comm paper spread';...
    '5yr rate'; '10yr rate'};
Titlenames = {'1yr rate';'CPI';'IP';'EBP'};
%% Plot the IRFs

idx = 1;
for ii=1:size(IRFs_single, 2)
    subplot(4,2,idx)
    xlim(subplot(4,2,idx),[min(min(IRFs_singleLo(:,ii)),...
        min(IRFSingle4Lo(:,ii))) max(max(IRFs_singleUp(:,ii), ...
        max(IRFSingle4Lo(:,ii))))]);
    plot(IRFs_single(:,ii),'-r','LineWidth',2); hold on; 
    plot(IRFs_singleLo(:,ii),'--r','LineWidth',1); hold on; 
    plot(IRFs_singleUp(:,ii),'--r','LineWidth',1); hold on; 
    plot(zeros(h),'-k')
    title(Titlenames{ii},'FontWeight','bold')
    axis tight
    idx = idx + 2;
end
legend('Original GK2015')
idx = 2;
for ii=1:4
    subplot(4,2,idx)
    xlim(subplot(4,2,idx),[min(min(IRFs_singleLo(:,ii)),...
        min(IRFSingle4Lo(:,ii))) max(max(IRFs_singleUp(:,ii), ...
        max(IRFSingle4Lo(:,ii))))]);
    plot(IRFSingle4(:,ii),'-k','LineWidth',2); hold on; 
    plot(IRFSingle4Up(:,ii),'--k','LineWidth',1); hold on; 
    plot(IRFSingle4Lo(:,ii),'--k','LineWidth',1); hold on; 
    plot(zeros(h),'-k')
    title(TitlenamesF{ii},'FontWeight','bold')
    axis tight
    idx = idx + 2;
end
legend('FAVAR with GK2015 + 4PC')
print('GK_FAV4_Comp','-dpdf','-r0');


idx = 1;
for ii=2:size(IRFSingle4, 2)
    subplot(3,2,idx) 
    plot(IRFSingle4(:,ii),'-k','LineWidth',2); hold on; 
    plot(IRFSingle4Up(:,ii),'--k','LineWidth',1); hold on; 
    plot(IRFSingle4Lo(:,ii),'--k','LineWidth',1); hold on; 
    plot(zeros(h),'-k')
    title(TitlenamesF{ii},'FontWeight','bold')
    axis tight
    idx = idx + 1;
end   
hold off
