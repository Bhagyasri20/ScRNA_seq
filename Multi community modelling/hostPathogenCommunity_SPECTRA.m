%% hostPathogenCommunity_SPECTRA_both.m
% MKN28 (host) + 7 AGORA2 microbes on Ham's F-12
%   shared       : all members share ONE extracellular pool [e]
%   compartment  : every member keeps its OWN [e]; a shared lumen [u]
clear; clc
initCobraToolbox(0); 
changeCobraSolver('gurobi','all'); 
addpath(genpath(pwd));

%% 1) settings
D1 = 'C:\Users\ommiv\Downloads\MKN28_models\';  D2 = 'C:\Users\ommiv\Downloads\files\';
files = [{[D1 'mkn28_0_updated.xml']}, strcat(D2,{ ...
  'Escherichia_coli_UTI89_UPEC_AGORA2_fixed.xml','Lactobacillus_kalixensis_DSM_16043_AGORA2_v2_01_fixed.xml', ...
  'Pseudomonas_aeruginosa_NCGM2_S1_AGORA2_fixed.xml','Prevotella_melaninogenica_ATCC_25845_AGORA2_fixed.xml', ...
  'Helicobacter_pylori_26695_AGORA2_fixed.xml','Streptococcus_salivarius_JIM8777_AGORA2_fixed.xml', ...
  'Fusobacterium_nucleatum_subsp_nucleatum_ATCC_25586_AGORA2_fixed.xml'})];
abbr  = {'MKN28','ecoli','lkal','paer','pmel','hpyl','ssal','fnuc'};      % first = host
% max uptake (mmol/gDW/h) = C(mmol/L)*V(mL) / (1000 * t(h) * DW(gDW/cell) * N(cells)), N*DW summed over ALL members
% No culture data: literature-typical ASSUMPTIONS, state them in Methods and test +/-2x.
Vol = 1;  t = 24;                          % mL of medium, hours
hostN = 5e5*Vol;   hostDW = 3e-10;         % MKN28 cells, ~300 pg dry weight/cell
MOI   = 10;        bacDW  = 3e-13;         % bacteria per host cell (per species), ~0.3 pg dry weight/bacterium
n     = numel(files);
totDW = hostN*hostDW + (n-1)*(MOI*hostN)*bacDW;   % total community biomass in Vol mL (gDW)
extra = {'pheme','sheme','mqn7','mqn8','2dmmq8','4abz','nac','orn_L','spmd','26dap_M','ura','xan','uri','cytd','adn','gsn','dgsn'};
trace = {'mn2','cobalt2','ni2','fe3','sel','nh4'};
gFrac = 0.5;  gCap = 0.1;                  % force each member to grow at >= min(gFrac*its own max, gCap) 1/h, so SPECTRA must keep a growth-capable network
hostExtra = {'orot','utp','dcyt','etha','ppa','acetone','cholate','ps_hs','pglyc_hs'};   % serum-like extras, added ONLY if the host cannot grow on F-12 alone
hostCap = Inf;                             % host growth cap (1/h). Inf = none, so time points can differ; a cap (e.g. 0.03) makes every time point equal
tpFiles = {'mkn28_0_updated','mkn28_2_updated','mkn28_6_updated'};  hrs = [0 2 6];   % host models, one per time point

% Ham's F-12, mM (Thermo Fisher 11765), salts as free ions, phenol red omitted
F12 = {'glc_D',10;'pyr',1;'gly',0.1;'ala_L',0.1;'arg_L',1;'asn_L',0.1;'asp_L',0.1;'cys_L',0.2;'glu_L',0.1;'gln_L',1; ...
 'his_L',0.1;'ile_L',0.031;'leu_L',0.1;'lys_L',0.2;'met_L',0.03;'phe_L',0.03;'pro_L',0.3;'ser_L',0.1;'thr_L',0.1; ...
 'trp_L',0.01;'tyr_L',0.03;'val_L',0.1;'btn',3e-5;'chol',0.1;'pnto_R',1e-3;'fol',2.9e-3;'ncam',2.95e-4;'pydxn',2.9e-4; ...
 'ribflv',9.8e-5;'thm',8.9e-4;'cbl1',1.03e-3;'inost',0.1;'ca2',0.3;'mg2',0.6;'k',2.98;'na1',148;'cl',137;'hco3',14; ...
 'pi',1;'fe2',3e-3;'zn2',3e-3;'cu2',1e-5;'so4',6e-3;'hxan',0.03;'thymd',2.9e-3;'ptrc',1e-3;'lipoate',1e-3;'lnlc',3e-4};
ions = {'na1','cl','k','ca2','mg2','pi','hco3','fe2','zn2','cu2','so4'};
vits = {'btn','pnto_R','fol','ncam','pydxn','ribflv','thm','cbl1'};

Mdl = cell(n,1);  bmId = cell(n,1);                                % load once, reuse for both builds
for i = 1:n
    Mdl{i} = readCbModel(files{i});
    bmId{i} = [abbr{i} '_' Mdl{i}.rxns{find(Mdl{i}.c,1)}];         % each model's own objective = biomass
end

modeList = 1:2;                                                  % 1 = shared pool, 2 = compartmentalised
useGapFill = true;  tkList = 1:numel(hrs);                     % switches: gap-fill lkal/pmel or not; which time points to run
gapAb  = {'lkal','pmel'};  if ~useGapFill, gapAb = {}; end                                        % members that cannot grow on their own -> gap-filled by SPECTRA
donors = find(~ismember(abbr,['MKN28',gapAb]));                  % the members that do grow: their reactions are the candidate pool
modes = {'shared','compartment'};  hostG = nan(numel(hrs),2);
for tk = tkList
Mdl{1} = readCbModel([D1 tpFiles{tk} '.xml']);                    % host model of this time point (same biomass reaction in all three)
for mi = modeList
compart = (mi == 2);  name = [modes{mi} '_' num2str(hrs(tk)) 'h'];  clear B
fprintf('\n================ %s model ================\n',name);

%% 2) prefix and merge (own exchanges dropped)
for i = 1:n
    m = Mdl{i};  ab = abbr{i};
    keep = ~startsWith(m.rxns,'EX_');  mets = m.mets(:);
    nE = ~endsWith(mets,'[e]') | compart;                          % compartment mode: [e] is prefixed too -> one [e] per member
    mets(nE) = strcat(ab,'_',mets(nE));
    B(i) = struct('S',m.S(:,keep),'mets',{mets},'rxns',{strcat(ab,'_',m.rxns(keep))},'lb',m.lb(keep),'ub',m.ub(keep)); %#ok<SAGROW>
end
% gap-fill candidates: every reaction of the growing members that lkal/pmel do not have yet, renamed for lkal/pmel
cand = {};  nB = n;
for g = 1:numel(gapAb)
    ti = find(strcmp(abbr,gapAb{g}));  taken = Mdl{ti}.rxns;
    for d = donors
        m = Mdl{d};
        keep = ~startsWith(m.rxns,{'EX_','DM_','sink_'}) & ~contains(lower(m.rxns),'biomass') & m.c == 0 & ~ismember(m.rxns,taken);
        mets = m.mets(:);  nE = ~endsWith(mets,'[e]') | compart;  mets(nE) = strcat(abbr{ti},'_',mets(nE));
        nB = nB + 1;
        B(nB) = struct('S',m.S(:,keep),'mets',{mets},'rxns',{strcat(abbr{ti},'_',m.rxns(keep))},'lb',min(m.lb(keep),0),'ub',m.ub(keep)); %#ok<SAGROW>
        taken = [taken; m.rxns(keep)];  cand = [cand; B(nB).rxns]; %#ok<AGROW>
    end
end
fprintf('gap-fill candidates: %d reactions\n',numel(cand));
[U,~,ic] = unique(vertcat(B.mets),'stable');
r0 = [0 cumsum(cellfun(@numel,{B.mets}))];  c0 = [0 cumsum(cellfun(@numel,{B.rxns}))];
I = []; J = []; V = [];
for i = 1:numel(B), [r,c,v] = find(B(i).S); I = [I; ic(r0(i)+r)]; J = [J; c0(i)+c]; V = [V; v]; end %#ok<AGROW>
Com.S = sparse(I,J,V,numel(U),c0(end));  Com.mets = U;
Com.rxns = vertcat(B.rxns);  Com.lb = vertcat(B.lb);  Com.ub = vertcat(B.ub);

if compart        % shared lumen [u] + one transporter  met[u] <=> member_met[e]  per member and metabolite
    ke = cell(n,1);  bs = cell(n,1);  lum = {};
    for i = 1:n
        ke{i} = find(startsWith(Com.mets,[abbr{i} '_']) & endsWith(Com.mets,'[e]') & full(sum(abs(Com.S),2)) > 0);
        bs{i} = strrep(extractAfter(Com.mets(ke{i}),[abbr{i} '_']),'[e]','');
        lum = [lum; bs{i}]; %#ok<AGROW>
    end
    lum = unique(lum);  nm0 = numel(Com.mets);
    Com.mets = [Com.mets; strcat(lum,'[u]')];  Com.S = [Com.S; sparse(numel(lum),size(Com.S,2))];
    for i = 1:n
        [~,li] = ismember(bs{i},lum);  k = numel(li);
        Com.S = [Com.S, sparse([nm0+li; ke{i}],[(1:k)'; (1:k)'],[-ones(k,1); ones(k,1)],numel(Com.mets),k)];
        Com.rxns = [Com.rxns; strcat(abbr{i},'_IEX_',bs{i},'[u]tr')];
        Com.lb = [Com.lb; -1000*ones(k,1)];  Com.ub = [Com.ub; 1000*ones(k,1)];
    end
    pool = '[u]';
else
    pool = '[e]';
end

e = find(endsWith(Com.mets,pool) & full(sum(abs(Com.S),2)) > 0);   % one exchange per pool metabolite = the medium
Com.S = [Com.S, sparse(e,1:numel(e),-1,numel(Com.mets),numel(e))];
Com.rxns = [Com.rxns; strcat('EX_',strrep(Com.mets(e),pool,''),'(',pool(2),')')];
Com.lb = [Com.lb; zeros(numel(e),1)];  Com.ub = [Com.ub; 1000*ones(numel(e),1)];
for q = find(contains(Com.mets,'biomass'))'                      % drain for biomass products nobody consumes
    if ~any(Com.S(q,:) < 0)
        Com.S(:,end+1) = sparse(q,1,-1,numel(Com.mets),1);  Com.rxns{end+1,1} = ['DM_' Com.mets{q}];
        Com.lb(end+1,1) = 0;  Com.ub(end+1,1) = 1000;
    end
end
Com.c = zeros(numel(Com.rxns),1);  Com.b = zeros(numel(Com.mets),1);  Com.csense = repmat('E',numel(Com.mets),1);
fprintf('community: %d reactions, %d metabolites, %d pool exchanges\n',numel(Com.rxns),numel(Com.mets),numel(e));

%% 3) medium: Ham's F-12 on the pool
for k = 1:size(F12,1)
    id = F12{k,1};
    if ismember(id,ions), v = 1000; elseif ismember(id,vits), v = 0.01; else, v = F12{k,2}*Vol/(1000*t*totDW); end
    Com = setLB(Com,id,v);
end
for id = {'h2o','h','co2'}, Com = setLB(Com,id{1},1000); end
Com = setLB(Com,'o2',5);
for id = trace, Com = setLB(Com,id{1},1); end                    % needed by AGORA2 biomass, not in F-12
for id = extra, Com = setLB(Com,id{1},1); end                    % heme / vitamin K etc. (state in Methods)
Com = setLB(Com,'crm_hs',0.01);                                  % serum-like ceramide: MKN28 cannot make sphingomyelin otherwise
if isfinite(hostCap), Com.ub(strcmp(Com.rxns,bmId{1})) = hostCap; end

%% 4) growth check on the community BEFORE SPECTRA, then force growth so it cannot be extracted away
bmIdx = zeros(n,1);  for i = 1:n, bmIdx(i) = find(strcmp(Com.rxns,bmId{i})); end
fba = maxGrowth(Com,bmIdx);
if fba(1) <= 1e-9                                                % host cannot grow on F-12 + supplements -> add serum-like extras
    for id = hostExtra, Com = setLB(Com,id{1},1); end
    fba = maxGrowth(Com,bmIdx);
    fprintf('host had no growth on F-12 alone; added extras: %s\n',strjoin(hostExtra,', '));
end
for i = 1:n, fprintf('%-6s growth in community BEFORE extraction = %.6g\n',abbr{i},fba(i)); end
grow = fba > 1e-9;
if ~all(grow), fprintf('NO growth possible (check medium / gap-fill) for: %s\n',strjoin(abbr(~grow),', ')); end
frac = gFrac;                                                    % joint feasibility: members compete for the same medium, so lower the forced rate until all fit
while true
    Z = Com;  Z.c(:) = 0;
    Z.lb(bmIdx(grow)) = min(frac*fba(grow),gCap);
    s = optimizeCbModel(Z,'max');
    if s.stat == 1, break; end
    frac = frac/2;
    if frac < 1e-4, error('Members cannot grow together on this medium; raise uptake bounds (smaller totDW) or add nutrients'); end
end
Com.lb(bmIdx(grow)) = Z.lb(bmIdx(grow));
fprintf('forced growth lower bounds (frac %.4g): %s\n',frac,strjoin(strcat(abbr(grow),'=',string(round(Com.lb(bmIdx(grow))',5))),', '));

%% 5) SPECTRA: consistency check + extraction in one step
core = unique([find(startsWith(Com.rxns,[abbr{1} '_'])); bmIdx(grow)]);   % whole host (incl. its transporters) + every growing member's biomass
core = setdiff(core,bmIdx(~grow));                               % a biomass that cannot carry flux must not be forced
w    = ones(numel(Com.rxns),1);  w(core) = 0;                    % non-core reactions weigh 1
isCand = ismember(Com.rxns,cand);  isTarget = false(size(w));
for g = 1:numel(gapAb), isTarget = isTarget | (startsWith(Com.rxns,[gapAb{g} '_']) & ~isCand); end
w(isTarget) = 0;                                                 % lkal/pmel's own reactions are free, added candidates cost 1 (as in the paper's gap-filling)
Final = spectraCCME(Com, core, 1e-5, 'stoichiometry', w, 1, {}, 'minNetLP', 7200);
if iscell(Final), Final = Final{1}; end

%% 6) check (lower bounds of biomass released, so this is the true maximum)
for i = 1:n
    j = find(strcmp(Final.rxns,bmId{i}));
    if isempty(j), fprintf('%-6s biomass removed (blocked on this medium)\n',abbr{i}); continue; end
    Y = Final;  Y.c(:) = 0;  Y.c(j) = 1;  Y.lb(j) = 0;  s = optimizeCbModel(Y,'max');
    fprintf('%-6s growth = %.6g  (solver status %d)\n',abbr{i},s.f,s.stat);
    if i == 1, hostG(tk,mi) = s.f; end
end
ok = spectraCC(Final,1e-5);
fprintf('%s: %d -> %d reactions, blocked after extraction: %d\n',name,numel(Com.rxns),numel(Final.rxns),numel(Final.rxns)-numel(ok));
added = intersect(cand,Final.rxns);                              % reactions SPECTRA had to borrow for lkal/pmel
for g = 1:numel(gapAb), fprintf('%s: %d reactions added from other members\n',gapAb{g},nnz(startsWith(added,[gapAb{g} '_']))); end
writetable(table(added),['added_reactions_' name '.csv']);
save(['Community_F12_' name '.mat'],'Com','Final','core','added','fba','-v7.3');
%% 7) export the extracted community model
CommunityModel = Final;                                          % the SPECTRA-extracted community of this run
save(['CommunityModel_' name '.mat'],'CommunityModel','-v7.3');  % MATLAB format, one file per model
try writeCbModel(Final,'format','sbml','fileName',['CommunityModel_' name '.xml']);   % SBML for other tools
catch ME, fprintf('SBML export skipped (%s)\n',ME.message); end
writetable(table(Final.rxns,printRxnFormula(Final,'printFlag',false),Final.lb,Final.ub,'VariableNames',{'rxn','formula','lb','ub'}),['CommunityModel_' name '_reactions.csv']);
fprintf('exported CommunityModel_%s (.mat, .xml, _reactions.csv)\n',name);
end
end
fprintf('\nMKN28 growth in the community (rows = %s h):\n',strjoin(string(hrs),'/'));
disp(array2table(hostG,'VariableNames',modes,'RowNames',strcat(string(hrs'),'h')));

function C = setLB(C,id,v)      % open uptake of one pool metabolite (skips ones not in the community)
    j = find(ismember(C.rxns,{['EX_' id '(e)'],['EX_' id '(u)']}));
    if ~isempty(j), C.lb(j) = -v; end
end
function g = maxGrowth(C,idx)  % max flux of each biomass reaction on its own (0 if infeasible)
    g = zeros(numel(idx),1);
    for q = 1:numel(idx)
        Y = C;  Y.c(:) = 0;  Y.c(idx(q)) = 1;  r = optimizeCbModel(Y,'max');
        if r.stat == 1 && isfinite(r.f), g(q) = r.f; end
    end
end