%% sanityCheckCommunity.m
clear; clc
initCobraToolbox(0);
changeCobraSolver('gurobi','all');
%% settings
names = {'shared_0h','compartment_0h','shared_2h','compartment_2h','shared_6h','compartment_6h'};
D1 = 'C:\Users\ommiv\Downloads\MKN28_models\';  D2 = 'C:\Users\ommiv\Downloads\files\';
files = [{[D1 'mkn28_0_updated.xml']}, strcat(D2,{ ...
  'Escherichia_coli_UTI89_UPEC_AGORA2_fixed.xml', ...
  'Lactobacillus_kalixensis_DSM_16043_AGORA2_v2_01_fixed.xml', ...
  'Pseudomonas_aeruginosa_NCGM2_S1_AGORA2_fixed.xml', ...
  'Prevotella_melaninogenica_ATCC_25845_AGORA2_fixed.xml', ...
  'Helicobacter_pylori_26695_AGORA2_fixed.xml', ...
  'Streptococcus_salivarius_JIM8777_AGORA2_fixed.xml', ...
  'Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586_AGORA2_fixed.xml'})];
abbr = {'MKN28','ecoli','lkal','paer','pmel','hpyl','ssal','fnuc'};
cfg.tol = 1e-6;            % flux above this = "carries flux"
cfg.fixCycles = true;      % remove energy-generating cycles
cfg.maxFix = 30;           % max blocking steps per member and cofactor
cfg.keepGrowth = 0.9;      % a fix is rejected if any member keeps < 90% of its growth
cfg.removeBlocked = true;  % delete blocked reactions at the end
cfg.checkStoich = true;    % stoichiometric consistency (slow, report only)


% expected medium = F-12 + the supplements
allowed = {'glc_D','pyr','gly','ala_L','arg_L','asn_L','asp_L','cys_L','glu_L','gln_L','his_L','ile_L','leu_L','lys_L','met_L','phe_L','pro_L', ...
  'ser_L','thr_L','trp_L','tyr_L','val_L','btn','chol','pnto_R','fol','ncam','pydxn','ribflv','thm','cbl1','inost','ca2','mg2','k','na1','cl', ...
  'hco3','pi','fe2','zn2','cu2','so4','hxan','thymd','ptrc','lipoate','lnlc','h2o','h','co2','o2','mn2','cobalt2','ni2','fe3','sel','nh4', ...
  'pheme','sheme','mqn7','mqn8','2dmmq8','4abz','nac','orn_L','spmd','26dap_M','ura','xan','uri','cytd','adn','gsn','dgsn','crm_hs', ...
  'orot','utp','dcyt','etha','ppa','acetone','cholate','ps_hs','pglyc_hs'};

% each test: a hydrolysis/oxidation demand that must carry NO flux when all uptakes are closed
T = {'atp',  {'atp[c]','h2o[c]','adp[c]','pi[c]','h[c]'},        [-1 -1 1 1 1];
     'gtp',  {'gtp[c]','h2o[c]','gdp[c]','pi[c]','h[c]'},        [-1 -1 1 1 1];
     'nadh', {'nadh[c]','nad[c]','h[c]'},                         [-1 1 1];
     'nadph',{'nadph[c]','nadp[c]','h[c]'},                       [-1 1 1];
     'fadh2',{'fadh2[c]','fad[c]','h[c]'},                        [-1 1 2];
     'accoa',{'accoa[c]','h2o[c]','coa[c]','ac[c]','h[c]'},       [-1 -1 1 1 1]};

% biomass reaction of every member = its own objective
bmId = cell(numel(files),1);
for i = 1:numel(files), m = readCbModel(files{i}); bmId{i} = [abbr{i} '_' m.rxns{find(m.c,1)}]; end
n = numel(abbr);  summary = {};

for q = 1:numel(names)
name = names{q};
fprintf('\n================ %s ================\n',name);
L = load(['CommunityModel_' name '.mat']);  M = L.CommunityModel;
if ~isfield(M,'b') || numel(M.b) ~= numel(M.mets), M.b = zeros(numel(M.mets),1); end
M.csense = repmat('E',numel(M.mets),1);  M.c = zeros(numel(M.rxns),1);
added = {};  if exist(['added_reactions_' name '.csv'],'file'), t = readtable(['added_reactions_' name '.csv']); added = t{:,1}; end
nStart = numel(M.rxns);  fixLog = {};

%% 1) integrity
nNaN = nnz(isnan(M.S));  badB = find(M.lb > M.ub);  infB = find(~isfinite(M.lb) | ~isfinite(M.ub));
dupR = numel(M.rxns) - numel(unique(M.rxns));  dupM = numel(M.mets) - numel(unique(M.mets));
emptyR = find(full(sum(abs(M.S),1)) == 0);  unusedM = find(full(sum(abs(M.S),2)) == 0);
fprintf('1) integrity: NaN in S %d | lb>ub %d | infinite bounds %d | duplicate rxn names %d | duplicate met names %d | empty rxns %d | unused mets %d\n', ...
    nNaN,numel(badB),numel(infB),dupR,dupM,numel(emptyR),numel(unusedM));
for j = badB(:)', [M.lb(j),M.ub(j)] = deal(min(M.lb(j),M.ub(j)),max(M.lb(j),M.ub(j))); end        % FIX: swap reversed bounds
M.lb(M.lb < -1000) = -1000;  M.ub(M.ub > 1000) = 1000;  M.lb(isinf(M.lb)) = -1000;  M.ub(isinf(M.ub)) = 1000;   % FIX: cap
if ~isempty(emptyR), M = removeRxns(M,M.rxns(emptyR)); end                                         % FIX: drop empty reactions
unusedM = find(full(sum(abs(M.S),2)) == 0);
if ~isempty(unusedM), M = removeMetabolites(M,M.mets(unusedM),false); end                          % FIX: drop unused metabolites
if dupR > 0 || dupM > 0, fprintf('   WARNING: duplicate names exist - fix them in the build script\n'); end

isEx = @(X) startsWith(X.rxns,'EX_');
bmIdx = @(X) cellfun(@(b) find(strcmp(X.rxns,b),1),bmId,'UniformOutput',false);

%% 2) open medium
ex = find(isEx(M));  openEx = ex(M.lb(ex) < 0);
ids = regexprep(M.rxns(openEx),'^EX_(.*)\([eu]\)$','$1');
extraOpen = openEx(~ismember(ids,allowed));
fprintf('2) medium: %d exchanges, %d with uptake allowed, %d not in F-12/supplement list\n',numel(ex),numel(openEx),numel(extraOpen));
if ~isempty(extraOpen), fprintf('   unexpected open uptakes: %s\n',strjoin(M.rxns(extraOpen)',', ')); end
writetable(table(M.rxns(openEx),M.lb(openEx),'VariableNames',{'exchange','lb'}),['sanity_open_medium_' name '.csv']);

%% 3) growth of each member
g0 = growthAll(M,bmId);
for i = 1:n, fprintf('3) growth %-6s = %.6g\n',abbr{i},g0(i)); end
if any(g0 <= 1e-9), fprintf('   WARNING: no growth for: %s\n',strjoin(abbr(g0 <= 1e-9),', ')); end

%% 4) free energy / cofactor generation with every uptake closed
nFlag = 0;  nFix = 0;
for k = 1:size(T,1)
  for i = 1:n
    ab = abbr{i};  it = 0;
    while true
        [f,x] = lunch(M,ab,T(k,:),cfg.tol);
        if isnan(f), break; end                                       % cofactor not present in this member
        if f <= cfg.tol, break; end
        if it == 0, nFlag = nFlag + 1; fprintf('4) %s makes %s from nothing: flux %.4g\n',ab,T{k,1},f); end
        if ~cfg.fixCycles || it >= cfg.maxFix, break; end
        e0 = isEx(M);  cyc = find(abs(x(1:numel(M.rxns))) > cfg.tol & ~e0);
        isAdd = ismember(M.rxns(cyc),added);  cyc = [cyc(isAdd); cyc(~isAdd)];      % blame gap-filled reactions first
        done = false;
        for c = cyc(:)'
            M2 = M;  if x(c) > 0, M2.ub(c) = 0; else, M2.lb(c) = 0; end              % block only the direction used by the cycle
            f2 = lunch(M2,ab,T(k,:),cfg.tol);
            if isnan(f2) || f2 >= f - 1e-9, continue; end
            g2 = growthAll(M2,bmId);
            if any(g0 > 0 & g2 < cfg.keepGrowth*g0), continue; end
            M = M2;  done = true;  nFix = nFix + 1;
            fixLog(end+1,:) = {name,ab,T{k,1},M.rxns{c},sign(x(c))}; %#ok<AGROW>
            break
        end
        if ~done, fprintf('   could not remove the %s/%s cycle automatically (%d reactions involved) - inspect sanity_cycle_%s.csv\n',ab,T{k,1},numel(cyc),name); ...
            writetable(table(M.rxns(cyc),x(cyc),'VariableNames',{'rxn','flux'}),['sanity_cycle_' name '.csv']); break; end
        it = it + 1;
    end
  end
end
fprintf('4) energy/cofactor tests: %d flagged, %d fixed by direction blocks\n',nFlag,nFix);
if ~isempty(fixLog), writetable(cell2table(fixLog,'VariableNames',{'model','member','test','blocked_rxn','direction'}),['sanity_cycle_fixes_' name '.csv']); end

%% 5) blocked reactions
ok = spectraCC(M,1e-5);  blocked = setdiff(1:numel(M.rxns),ok);
fprintf('5) blocked reactions: %d of %d\n',numel(blocked),numel(M.rxns));
for i = 1:n, fprintf('   %-6s blocked %d\n',abbr{i},nnz(startsWith(M.rxns(blocked),[abbr{i} '_']))); end
if cfg.removeBlocked && ~isempty(blocked)
    M2 = removeRxns(M,M.rxns(blocked));
    un = find(full(sum(abs(M2.S),2)) == 0);  if ~isempty(un), M2 = removeMetabolites(M2,M2.mets(un),false); end
    g2 = growthAll(M2,bmId);
    if any(g0 > 0 & g2 < cfg.keepGrowth*g0), fprintf('   removal would cut growth - blocked reactions kept\n');
    else, M = M2;  fprintf('   FIX: removed %d blocked reactions\n',numel(blocked)); end
end

%% 6) dead-end metabolites
Sp = M.S > 0;  Sn = M.S < 0;
canProd = full(sum(Sp(:, M.ub > 0),2) > 0) | full(sum(Sn(:, M.lb < 0),2) > 0);
canCons = full(sum(Sn(:, M.ub > 0),2) > 0) | full(sum(Sp(:, M.lb < 0),2) > 0);
dead = find(~(canProd & canCons));
fprintf('6) dead-end metabolites: %d of %d\n',numel(dead),numel(M.mets));
writetable(table(M.mets(dead),'VariableNames',{'metabolite'}),['sanity_deadends_' name '.csv']);

%% 7) stoichiometric consistency
nInc = NaN;
if cfg.checkStoich
    try
        M.SIntRxnBool = ~isEx(M) & ~startsWith(M.rxns,'DM_') & ~contains(M.rxns,'_IEX_');
        [~,~,incM] = findStoichConsistentSubset(M,0,0);  nInc = nnz(incM);
        fprintf('7) stoichiometrically inconsistent metabolites: %d\n',nInc);
    catch ME, fprintf('7) stoichiometry check skipped (%s)\n',ME.message); end
end

%% 8) community solution: mass balance + cross-feeding
Y = M;  Y.c(:) = 0;
for i = 1:n, j = find(strcmp(Y.rxns,bmId{i}),1);  if ~isempty(j), Y.c(j) = 1; end, end
s = optimizeCbModel(Y,'max','one');
if s.stat == 1
    v = s.x;  res = max(abs(M.S*v));
    fprintf('8) community max total growth %.6g | max |S*v| = %.3g\n',s.f,res);
    pool = find(endsWith(M.mets,{'[e]','[u]'}));  rows = {};
    for r = pool(:)'
        net = zeros(n,1);
        for i = 1:n
            cols = find(startsWith(M.rxns,[abbr{i} '_']) & M.S(r,:)' ~= 0);
            net(i) = sum(full(M.S(r,cols)') .* v(cols));
        end
        pr = find(net > cfg.tol);  co = find(net < -cfg.tol);
        for a = pr(:)', for b = co(:)'
            rows(end+1,:) = {M.mets{r},abbr{a},abbr{b},min(net(a),-net(b))}; %#ok<AGROW>
        end, end
    end
    if ~isempty(rows), writetable(cell2table(rows,'VariableNames',{'metabolite','producer','consumer','flux'}),['sanity_crossfeeding_' name '.csv']); end
    fprintf('   cross-fed metabolites (producer -> consumer pairs): %d  (sanity_crossfeeding_%s.csv)\n',size(rows,1),name);
else
    fprintf('8) community LP not solved (status %d)\n',s.stat);
end

%% save the checked model
g1 = growthAll(M,bmId);  CommunityModelChecked = M;
save(['CommunityModel_' name '_checked.mat'],'CommunityModelChecked','-v7.3');
try writeCbModel(M,'format','sbml','fileName',['CommunityModel_' name '_checked.xml']); catch, end
fprintf('final: %d -> %d reactions; growth after fixes: %s\n',nStart,numel(M.rxns),strjoin(strcat(abbr,'=',string(round(g1',4))),', '));
summary(end+1,:) = {name,nStart,numel(M.rxns),nFlag,nFix,numel(blocked),numel(dead),nInc,min(g1)}; %#ok<SAGROW>
end
S = cell2table(summary,'VariableNames',{'model','rxns_before','rxns_after','energy_flags','energy_fixed','blocked','deadend_mets','inconsistent_mets','min_growth'});
disp(S);  writetable(S,'SanityReport.csv');

%% helpers
function g = growthAll(M,bmId)       % max growth of each member on the model's own medium (NaN if its biomass is absent)
    g = nan(numel(bmId),1);
    for i = 1:numel(bmId)
        j = find(strcmp(M.rxns,bmId{i}),1);  if isempty(j), continue; end
        Y = M;  Y.c(:) = 0;  Y.c(j) = 1;  Y.lb(j) = 0;  s = optimizeCbModel(Y,'max');
        g(i) = 0;  if s.stat == 1 && isfinite(s.f), g(i) = s.f; end
    end
end
function [f,x] = lunch(M,ab,t,tol)   % flux of a demand (e.g. ATP hydrolysis) with ALL uptakes closed; NaN if metabolites missing
    f = NaN;  x = [];
    mets = strcat(ab,'_',t{2});
    if ~all(ismember(mets,M.mets)), return; end
    Y = M;  ex = startsWith(Y.rxns,'EX_');  Y.lb(ex) = 0;  nR = numel(Y.rxns);
    Y = addReaction(Y,'TEST_DEMAND','metaboliteList',mets,'stoichCoeffList',t{3},'reversible',false,'lowerBound',0,'upperBound',1000);
    Y.c = zeros(numel(Y.rxns),1);  Y.c(end) = 1;
    s = optimizeCbModel(Y,'max');
    f = 0;  if s.stat == 1 && isfinite(s.f), f = s.f; end
    if nargout > 1 && f > tol, s = optimizeCbModel(Y,'max','one');  if s.stat == 1, x = s.x(1:nR); end, end
end
