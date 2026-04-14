%% SECTION 1: load genome scale yeast model
load yeast-GEM.mat;

%check S matrix (metabolite, reaction)
model.S(1,:);

%check reaction and metabolite names
model.rxnNames(3413);
model.metNames(1);
find(strcmp(model.metNames,'alpha-D-Glucose'));

%% SECTION 2: Let's find branch point metabolites: those that participate in 2 or more
%% than 2 consuming reactions, not counting exchange and growth reactions.
Nmets = size(model.S,1);
branchMets = [];
for i=1:Nmets
    myRxns = find(model.S(i,:));
    myConsumingRxns = find(model.S(i,:)<0);
    myRxnsNonExchange = [];
    nRxns = size(myRxns,2);
    nConsumptionRxns = size(myConsumingRxns,2);
    for j=1:nConsumptionRxns;
        if (strcmp(model.subSystems{myConsumingRxns(j)},'Exchange reaction') == false & strcmp(model.subSystems{myConsumingRxns(j)},'Growth') == false)
            myRxnsNonExchange = [myRxnsNonExchange,myConsumingRxns(j)];
        end
    end
    sumStoichs = sum(model.S(i,myRxnsNonExchange),2);
    if mod(sumStoichs,1) ~= 0;
        disp([num2str(i),'is decimal']);
    else 
        disp([num2str(i),'is integer']);
        if sumStoichs <= -2 & strcmp(model.metNames{i},'acetyl-CoA') == false ...
                & strcmp(model.metNames{i},'AMP') == false & strcmp(model.metNames{i},'ATP') == false ...
                & strcmp(model.metNames{i},'ADP') == false & strcmp(model.metNames{i},'NADH') == false ...
                & strcmp(model.metNames{i},'GTP') == false ...
                & strcmp(model.metNames{i},'GMP') == false ...
                & strcmp(model.metNames{i},'GDP') == false ...
                & strcmp(model.metNames{i},'FMN') == false ...
                & strcmp(model.metNames{i},'FMNH2') == false ...
                & strcmp(model.metNames{i},'FAD') == false ...
                & strcmp(model.metNames{i},'FADH2') == false ...
                & strcmp(model.metNames{i},'NAD') == false ...
                & strcmp(model.metNames{i},'NADP') == false ...
                & strcmp(model.metNames{i},'NADP(+)') == false ...
                & strcmp(model.metNames{i},'NADPH') == false ...
                & strcmp(model.metNames{i},'coenzyme A') == false ...
                & strcmp(model.metNames{i},'polyphosphate') == false ...
                & strcmp(model.metNames{i},'phosphate') == false ...
                & strcmp(model.metNames{i},'diphosphate') == false ...
                & strcmp(model.metNames{i},'ferricytochrome c') == false ...
                & strcmp(model.metNames{i},'ferrocytochrome c') == false ...
                & strcmp(model.metNames{i},'Ferricytochrome b5') == false ...
                & strcmp(model.metNames{i},'Ferrocytochrome b5') == false ...
                & strcmp(model.metNames{i},'iron(3+)') == false ...
                & strcmp(model.metNames{i},'iron(2+)') == false ...
                & strcmp(model.metNames{i},'superoxide') == false ...
                & strcmp(model.metNames{i},'hydrogen peroxide') == false ...
                & strcmp(model.metNames{i},'Ca(2+)') == false ...
                & strcmp(model.metNames{i},'glutathione') == false ...
                & strcmp(model.metNames{i},'quinone') == false ...
                & strcmp(model.metNames{i},'ubiquinone-6') == false ...
                & strcmp(model.metNames{i},'H2O') == false ...
                & strcmp(model.metNames{i},'H+') == false ...
                & strcmp(model.metNames{i},'oxygen') == false ...
                & strcmp(model.metNames{i},'carbon dioxide') == false ...
                & strcmp(model.metNames{i},'bicarbonate') == false ...
                & strcmp(model.metNames{i},'sodium') == false ...
                & strcmp(model.metNames{i},'ammonium') == false ...
                & strcmp(model.metNames{i},'TRX1') == false ...
                & strcmp(model.metNames{i},'TRX1 disulphide') == false ...;
            branchMets = [branchMets,i];
        end
    end
end
nBranchMets = size(branchMets,2);
branchMetNames = {};
branchMetNamesCytosol = {};
branchMetsCytosol = [];
branchMetNamesMito = {};
branchMetsMito = [];
branchMetCompsCounts = zeros(1,size(model.compNames,1));
for dum=1:nBranchMets;
    C = {};
    C={model.metNames{branchMets(dum)},num2str(model.metComps(branchMets(dum)))};
    newName = strjoin(C,'_');
    branchMetNames = [branchMetNames;newName];
    if model.metComps(branchMets(dum)) == 1;
        branchMetNamesCytosol = [branchMetNamesCytosol;newName];
        branchMetsCytosol = [branchMetsCytosol,branchMets(dum)];
    end
    if model.metComps(branchMets(dum)) == 9;
        branchMetNamesMito = [branchMetNamesMito;newName];
        branchMetsMito = [branchMetsMito,branchMets(dum)];
    end
    branchMetCompsCounts(model.metComps(branchMets(dum)))=branchMetCompsCounts(model.metComps(branchMets(dum)))+1;
end

%% SECTION 3: now that we got branch point metabolites, let's see which ones involve
%% NADH or ATP at the branch
branchMetsCoSubNADH = [];
branchMetsCoSubATP = [];
branchMetsCoSubNADPH = [];
for i=1:size(branchMets,2);
    myRxns = find(model.S(branchMets(i),:));
    myConsumingRxns = find(model.S(branchMets(i),:)<0);
    nConsumptionRxns = size(myConsumingRxns,2);
    dumb = 0;
    dumber = 0;
    dumberer=0;
    for j=1:nConsumptionRxns;
        consumedSubstrates = find(model.S(:,myConsumingRxns(j))<0)
        nSubstrates = size(consumedSubstrates,1);
        for jj=1:nSubstrates;
            if (strcmp(model.metNames{consumedSubstrates(jj)},'NAD') == true | strcmp(model.metNames{consumedSubstrates(jj)},'NADH') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NADH ']);
                dumb = 1;
            end
            if (strcmp(model.metNames{consumedSubstrates(jj)},'NADP(+)') == true | strcmp(model.metNames{consumedSubstrates(jj)},'NADPH') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NAD(P)H ']);
                dumberer = 1;
            end
            if (strcmp(model.metNames{consumedSubstrates(jj)},'ATP') == true | strcmp(model.metNames{consumedSubstrates(jj)},'ADP') == true | strcmp(model.metNames{consumedSubstrates(jj)},'AMP') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses ATP ']);
                dumber = 1;
            end
        end
    end
    if dumb == 1;
        branchMetsCoSubNADH = [branchMetsCoSubNADH,branchMets(i)];
    end
    if dumber == 1;
        branchMetsCoSubATP = [branchMetsCoSubATP,branchMets(i)];
    end
    if dumberer == 1;
        branchMetsCoSubNADPH = [branchMetsCoSubNADPH,branchMets(i)];
    end
end
tots=[size(branchMetsCoSubNADH),size(branchMetsCoSubATP),size(branchMetsCoSubNADPH)]

%% SECTION 4: let's do the same, but for cytosol metabolites only
%% but also chasing any cytosol metabolites going into mitochon.
branchMetsCoSubNADH_cytosol = [];
branchMetsCoSubATP_cytosol = [];
branchMetsCoSubNADPH_cytosol = [];
branchMetsCoSubMULTI_cytosol = [];
for i=1:size(branchMetsCytosol,2);
    myRxns = find(model.S(branchMetsCytosol(i),:)); 
    myConsumingRxns = find(model.S(branchMetsCytosol(i),:)<0);
    considerMito = 0;
    nConsumptionRxns = size(myConsumingRxns,2);
    for ii=1:nConsumptionRxns;
        if (strcmp(model.subSystems{myConsumingRxns(ii)},'Transport [c, m]') == true);
            considerMito = 1;
        end
    end
    if (considerMito == 1);
        metaboliteName = model.metNames{branchMetsCytosol(i)};
        metaboliteIDs = find(strcmp(model.metNames,metaboliteName));
        metaboliteComps = model.metComps(metaboliteIDs);
        metaboliteIndex = find(metaboliteComps==9);
        metaboliteMitID =  metaboliteIDs(metaboliteIndex);
        myRxns2 = find(model.S(metaboliteMitID,:)); 
        myConsumingRxns2 = find(model.S(metaboliteMitID,:)<0);
        myConsumingRxns=[myConsumingRxns,myConsumingRxns2];
        nConsumptionRxns = size(myConsumingRxns,2);
    end
    dumb = 0;
    dumber = 0;
    dumberer=0;

    for j=1:nConsumptionRxns;
        consumedSubstrates = find(model.S(:,myConsumingRxns(j))<0)
        nSubstrates = size(consumedSubstrates,1);
        for jj=1:nSubstrates;
            if (strcmp(model.metNames{consumedSubstrates(jj)},'NAD') == true | strcmp(model.metNames{consumedSubstrates(jj)},'NADH') == true)
                disp(['met ',num2str(branchMetsCytosol(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NADH ']);
                dumb = 1;

            end
            if (strcmp(model.metNames{consumedSubstrates(jj)},'NADP(+)') == true | strcmp(model.metNames{consumedSubstrates(jj)},'NADPH') == true)
                disp(['met ',num2str(branchMetsCytosol(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NAD(P)H ']);
                dumberer = 1;
            end
            if (strcmp(model.metNames{consumedSubstrates(jj)},'ATP') == true | strcmp(model.metNames{consumedSubstrates(jj)},'ADP') == true | strcmp(model.metNames{consumedSubstrates(jj)},'AMP') == true)
                disp(['met ',num2str(branchMetsCytosol(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses ATP ']);
                dumber = 1;
            end
        end
    end
    if dumb == 1 && dumber ==0 && dumberer==0;
        branchMetsCoSubNADH_cytosol = [branchMetsCoSubNADH_cytosol,branchMetsCytosol(i)];
    end
    if dumb == 0 && dumber ==1 && dumberer==0;
        branchMetsCoSubATP_cytosol = [branchMetsCoSubATP_cytosol,branchMetsCytosol(i)];
    end
    if dumb == 0 && dumber ==0 && dumberer==1;
        branchMetsCoSubNADPH_cytosol = [branchMetsCoSubNADPH_cytosol,branchMetsCytosol(i)];
    end
    if (dumb ==1 && dumber==1 && dumberer==1) || (dumb==1 && dumber==1 && dumberer==0 ) || (dumb==1 && dumber==0 && dumberer==1) || (dumb==0 && dumber==1 && dumberer==1) 
        branchMetsCoSubMULTI_cytosol=[branchMetsCoSubMULTI_cytosol,branchMetsCytosol(i)];
    end
end
model.metNames{branchMetsCoSubNADH_cytosol}
tots_cyto=[size(branchMetsCoSubNADH_cytosol),size(branchMetsCoSubATP_cytosol),size(branchMetsCoSubNADPH_cytosol),size(branchMetsCoSubMULTI_cytosol)]
cyto_metab_tot=sum(model.metComps(:)==1)
tots;
tots_cyto
size(branchMetsCytosol)
size(unique([branchMetsCoSubNADH_cytosol,branchMetsCoSubATP_cytosol,branchMetsCoSubNADPH_cytosol,branchMetsCoSubMULTI_cytosol]))


