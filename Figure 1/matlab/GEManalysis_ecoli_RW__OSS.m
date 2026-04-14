%% SECTION1: load genome scale E.coli model (downloaded from BIGG database: http://bigg.ucsd.edu/models/iML1515)
load iML1515.mat;

%check S matrix (metabolite, reaction)
iML1515.S(1,:);

%check reaction and metabolite names
iML1515.metNames(:)
iML1515.mets(:)
iML1515.rxnNames(:)

%% SECTION2: just checking names
for dum=1:1877
    if (contains(iML1515.metNames{dum},'Nicotinamide adenine dinucleotide phosphate')==true)
        display(dum);
    end
end

% SECTION3: Let's find branch point metabolites: those that participate in 2 or more
% than 2 consuming reactions, not counting exchange reactions.
Nmets = size(iML1515.S,1)
branchMets = [];
for i=1:Nmets
    myRxns = find(iML1515.S(i,:));
    myConsumingRxns = find(iML1515.S(i,:)<0);
    myRxnsNonExchange = [];
    nRxns = size(myRxns,2);
    nConsumptionRxns = size(myConsumingRxns,2);
    for j=1:nConsumptionRxns;
        if (strcmp(iML1515.subSystems{myConsumingRxns(j)},'Exchange reaction') == false)
            myRxnsNonExchange = [myRxnsNonExchange,myConsumingRxns(j)];
        end
    end
    sumStoichs = sum(iML1515.S(i,myRxnsNonExchange),2);
    if mod(sumStoichs,1) ~= 0;
        disp([num2str(i),'is decimal'])
    else 
        disp([num2str(i),'is integer'])
        if sumStoichs <= -2 & strcmp(iML1515.metNames{i},'Acetyl-CoA') == false ...
                & strcmp(iML1515.metNames{i},'AMP C10H12N5O7P') == false & strcmp(iML1515.metNames{i},'ATP C10H12N5O13P3TP') == false ...
                & strcmp(iML1515.metNames{i},'ADP C10H12N5O10P2P') == false & strcmp(iML1515.metNames{i},'Nicotinamide adenine dinucleotide - reduced') == false ...
                & strcmp(iML1515.metNames{i},'ADP C10H12N5O10P2') == false ...
                & strcmp(iML1515.metNames{i},'Nicotinamide adenine dinucleotide') == false ...
                & strcmp(iML1515.metNames{i},'Nicotinamide adenine dinucleotide phosphate') == false ...
                & strcmp(iML1515.metNames{i},'Nicotinamide adenine dinucleotide phosphate - reduced') == false ...
                & strcmp(iML1515.metNames{i},'Coenzyme A') == false ...
                & strcmp(iML1515.metNames{i},'Hydrogen peroxide') == false ...
                & strcmp(iML1515.metNames{i},'Chloride') == false ...
                & strcmp(iML1515.metNames{i},'Cu+') == false ...
                & strcmp(iML1515.metNames{i},'H2O H2O') == false ...
                & strcmp(iML1515.metNames{i},'Pyridoxine') == false ...
                & strcmp(iML1515.metNames{i},'Phosphate') == false ...
                & strcmp(iML1515.metNames{i},'CO2 CO2') == false ...
                & strcmp(iML1515.metNames{i},'Sodium') == false ...
                & strcmp(iML1515.metNames{i},'Fe(III)hydroxamate') == false ...
                & strcmp(iML1515.metNames{i},'Fe(III)hydoxamate, unloaded') == false ...
                & strcmp(iML1515.metNames{i},'Potassium') == false ...
                & strcmp(iML1515.metNames{i},'Nitric oxide') == false ...
                & strcmp(iML1515.metNames{i},'Iron (Fe3+)') == false ...
                & strcmp(iML1515.metNames{i},'Hydrogen') == false ...
                & strcmp(iML1515.metNames{i},'Reduced glutathione') == false ...
                & strcmp(iML1515.metNames{i},'Ferrichrome minus Fe(III)') == false ...
                & strcmp(iML1515.metNames{i},'IscU scaffold protein') == false ...
                & strcmp(iML1515.metNames{i},'IscS with bound sulfur') == false ...
                & strcmp(iML1515.metNames{i},'IscU with bound [2Fe-2S] cluster') == false ...
                & strcmp(iML1515.metNames{i},'Reduced riboflavin') == false ...
                & strcmp(iML1515.metNames{i},'Reduced thioredoxin') == false ...
                & strcmp(iML1515.metNames{i},'MoaD Protein with thiocarboxylate') == false ...
                & strcmp(iML1515.metNames{i},'Superoxide anion') == false ...
                & strcmp(iML1515.metNames{i},'Nickel') == false ...
                & strcmp(iML1515.metNames{i},'SufBCD scaffold complex') == false ...
                & strcmp(iML1515.metNames{i},'Oxidized glutathione') == false ...
                & strcmp(iML1515.metNames{i},'Ubiquinone-8') == false ...
                & strcmp(iML1515.metNames{i},'Ferrichrome') == false ...
                & strcmp(iML1515.metNames{i},'O2 O2') == false ...
                & strcmp(iML1515.metNames{i},'Ammonium') == false ...
                & strcmp(iML1515.metNames{i},'Reduced FMN') == false ...
                & strcmp(iML1515.metNames{i},'FMN C17H19N4O9P') == false ...
                & strcmp(iML1515.metNames{i},'Hydrogen peroxide') == false ...
                & strcmp(iML1515.metNames{i},'Flavodoxin semi oxidized') == false ...
                & strcmp(iML1515.metNames{i},'Flavodoxin reduced') == false ...
                & strcmp(iML1515.metNames{i},'Magnesium') == false ...
                & strcmp(iML1515.metNames{i},'Bicarbonate') == false ...
                & strcmp(iML1515.metNames{i},'Manganese') == false ...
                & strcmp(iML1515.metNames{i},'Periplasmic disulfide isomerase/thiol-disulphide oxidase (oxidized)') == false ...
                & strcmp(iML1515.metNames{i},'Periplasmic protein disulfide isomerase I (reduced)') == false ...
                & strcmp(iML1515.metNames{i},'Generic ferrioxamine-Fe-III') == false ...
                & strcmp(iML1515.metNames{i},'Fused thiol:disulfide interchange protein (reduced)') == false ...
                & strcmp(iML1515.metNames{i},'Protein disulfide isomerase II (oxidized)') == false ...
                & strcmp(iML1515.metNames{i},'Mercury  charged 2  Hg') == false ...
                & strcmp(iML1515.metNames{i},'Generic ferrioxamine-Fe-III') == false ...
                & strcmp(iML1515.metNames{i},'H+') == false ...
                & strcmp(iML1515.metNames{i},'Zinc') == false ...
                & strcmp(iML1515.metNames{i},'Cesium ion') == false ...
                & strcmp(iML1515.metNames{i},'Fe2+ mitochondria') == false ...
                & strcmp(iML1515.metNames{i},'Cadmium') == false ...
                & strcmp(iML1515.metNames{i},'Co2+') == false ...
                & strcmp(iML1515.metNames{i},'Copper') == false ...
                & strcmp(iML1515.metNames{i},'SufBCD with bound [2Fe-2S] cluster') == false ...
                & strcmp(iML1515.metNames{i},'Calcium') == false ...
                & strcmp(iML1515.metNames{i},'Coprogen unloaded (no Fe(III))') == false ...
                & strcmp(iML1515.metNames{i},'Coprogen') == false ...
                & strcmp(iML1515.metNames{i},'Ferroxamine minus Fe(3)') == false ...
                & strcmp(iML1515.metNames{i},'Glutaredoxin (reduced)') == false ...
                & strcmp(iML1515.metNames{i},'Flavin adenine dinucleotide reduced') == false ...;
            branchMets = [branchMets,i];
        end
    end
end
branchMets
nBranchMets = size(branchMets,2)
branchMetNames = {};
branchMetNamesCytosol = {};
branchMetsCytosol = [];
branchMetNamesPeriplasm = {};
branchMetsPeriplasm = [];
for dum=1:nBranchMets;
    branchMetNames = [branchMetNames;iML1515.metNames{branchMets(dum)}];
    if contains(iML1515.mets(branchMets(dum)),'_c') == true;
        branchMetNamesCytosol = [branchMetNamesCytosol;iML1515.metNames{branchMets(dum)}];
        branchMetsCytosol = [branchMetsCytosol,branchMets(dum)];
    end
    if contains(iML1515.mets(branchMets(dum)),'_p') == true;
        branchMetNamesPeriplasm = [branchMetNamesPeriplasm;iML1515.metNames{branchMets(dum)}];
        branchMetsPeriplasm = [branchMetsPeriplasm,branchMets(dum)];
    end
end

%% SECTION4: now that we got branch point metabolites, let's see which ones involve
%% NADH or ATP at the branch
branchMetsCoSubNADH = [];
branchMetsCoSubATP = [];
branchMetsCoSubNADPH = [];
for i=1:size(branchMets,2)
    myRxns = find(iML1515.S(branchMets(i),:));
    myConsumingRxns = find(iML1515.S(branchMets(i),:)<0)
    nConsumptionRxns = size(myConsumingRxns,2);
    dumb = 0;
    dumber = 0;
    dumberer = 0;
    for j=1:nConsumptionRxns;
        consumedSubstrates = find(iML1515.S(:,myConsumingRxns(j))<0)
        nSubstrates = size(consumedSubstrates,1);
        for jj=1:nSubstrates;
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide - reduced') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NADH ']);
                dumb = 1;
            end
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide phosphate') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide phosphate - reduced') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NAD(P)H ']);
                dumberer = 1;
            end
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'ATP C10H12N5O13P3') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'ADP C10H12N5O10P2') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'AMP C10H12N5O7P') == true)
                disp(['met ',num2str(branchMets(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses ATP ']);
                dumber = 1;
            end
        end
    end
    if dumb == 1 ;
        branchMetsCoSubNADH = [branchMetsCoSubNADH,branchMets(i)];
    end
    if dumber == 1;
        branchMetsCoSubATP = [branchMetsCoSubATP,branchMets(i)];
    end
    if dumberer == 1;
        branchMetsCoSubNADPH = [branchMetsCoSubNADPH,branchMets(i)];
    end
end
tots=[size(branchMetsCoSubNADH),size(branchMetsCoSubATP),size(branchMetsCoSubNADPH)];


%% SECTION5: let's do the same for cytosolic branch metabolites
branchMetsCoSubNADH_cytosol = [];
branchMetsCoSubATP_cytosol = [];
branchMetsCoSubNADPH_cytosol = [];
branchMetsCoSubMULTI_cytosol = [];
for i=1:size(branchMetsCytosol,2)
    myRxns = find(iML1515.S(branchMetsCytosol(i),:));
    myConsumingRxns = find(iML1515.S(branchMetsCytosol(i),:)<0)
    nConsumptionRxns = size(myConsumingRxns,2);
    dumb = 0;
    dumber = 0;
    dumberer = 0;
    for j=1:nConsumptionRxns;
        consumedSubstrates = find(iML1515.S(:,myConsumingRxns(j))<0)
        nSubstrates = size(consumedSubstrates,1);
        for jj=1:nSubstrates;
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide - reduced') == true)
                disp(['met ',num2str(branchMetsCytosol(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NADH ']);
                dumb = 1;
            end
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide phosphate') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'Nicotinamide adenine dinucleotide phosphate - reduced') == true)
                disp(['met ',num2str(branchMetsCytosol(i)),' rxn ',num2str(myConsumingRxns(j)),' consumption uses NAD(P)H ']);
                dumberer = 1;
            end
            if (strcmp(iML1515.metNames{consumedSubstrates(jj)},'ATP C10H12N5O13P3') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'ADP C10H12N5O10P2') == true | strcmp(iML1515.metNames{consumedSubstrates(jj)},'AMP C10H12N5O7P') == true)
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
iML1515.metNames{branchMetsCoSubNADH_cytosol}

tots_cyto=[size(branchMetsCoSubNADH_cytosol), size(branchMetsCoSubATP_cytosol), size(branchMetsCoSubNADPH_cytosol),size(branchMetsCoSubMULTI_cytosol)];

tots;
tots_cyto

tot_cyto_metabs=sum(contains(iML1515.mets(:),'_c') == true)
size(branchMetsCytosol)


size(unique([branchMetsCoSubNADH_cytosol,branchMetsCoSubATP_cytosol,branchMetsCoSubNADPH_cytosol]))


