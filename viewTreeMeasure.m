function viewTreeMeasure()

colors = distinguishable_colors(20);
colors([3,4,5,6,7,8,9,10,11],:) = []; 
colors = [[0.8, 0.8,0.8] ; colors];
set(0,'DefaultAxesColorOrder',colors);
set(0,'DefaultAxesLineStyleOrder',{'-',':','--','.-'});

%load('treeResults_with_newKang100.mat');
%load('fullTreeResultsTriplets.mat');
 load('simpleMeasurementTree.mat');

cahoy_human = load('cahoy_human_homolgues.mat');
[cahoy_human_subset.scores_neurons, cahoy_human_subset.neurons_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.neurons_human_entrez, cahoy_human.neurons_human_hsbc);
[cahoy_human_subset.scores_astro, cahoy_human_subset.astro_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.astro_human_entrez, cahoy_human.astro_human_hsbc);
[cahoy_human_subset.scores_oligo, cahoy_human_subset.oligo_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.oligo_human_entrez, cahoy_human.oligo_human_hsbc);

[housekeeping.ensembl_gene_id,housekeeping.transcript_gene_id ,housekeeping.entrez ,housekeeping.symbol ] = textread('/cortex/data/gene_sets/house_keeping_genes/Erez_2003.txt','%s %s %d %s','headerlines',1);
[housekeeping_subset.scores, housekeeping_subset.symbols, housekeeping_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, housekeeping.entrez, housekeeping.symbol);

[axon_guidance.symbol ] = textread('/cortex/data/gene_sets/human_axon_gudiance','%s');
[axon_guidance_subset.scores, axon_guidance_subset.symbols, axon_guidance_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, nan(size(axon_guidance.symbol)), axon_guidance.symbol);

figure;
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
set(scatterDots,'SizeData',20);
xlabel('Human6 triplets scores', 'fontsize',20); ylabel('Kang triplets scores', 'fontsize',20);


hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.neurons_human_entrez, cahoy_human.neurons_human_hsbc, colors(2,:));
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.astro_human_entrez, cahoy_human.astro_human_hsbc, colors(3,:));
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.oligo_human_entrez, cahoy_human.oligo_human_hsbc, colors(4,:));
%scatterScores(human6RandomResults, human_gene_info, brainspanRandomResults , developing_genes_info);

hleg = legend('All','Neurons','Astrocytes','Oligodendrocytes','Random');  set(hleg,'Location','SouthEast'); legend('boxoff');     
set(gca,'box','off');  
set(gca,'ytick',[0:0.1:1]); set(gca,'xtick',[0:0.1:1]);
title('');

figure;
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
set(scatterDots,'SizeData',20);
xlabel('Human6 triplets scores', 'fontsize',20); ylabel('Kang triplets scores', 'fontsize',20);
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, housekeeping.entrez ,housekeeping.symbol, [1.0,0.4,0.0] );
hleg = legend('All','House keeping genes');  set(hleg,'Location','SouthEast'); legend('boxoff');     
set(gca,'box','off');  
set(gca,'ytick',[0:0.1:1]); set(gca,'xtick',[0:0.1:1]);
title('');


figure;
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
set(scatterDots,'SizeData',20);
xlabel('Human6 triplets scores', 'fontsize',20); ylabel('Kang triplets scores', 'fontsize',20);
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan(size(axon_guidance.symbol )) ,axon_guidance.symbol , [1.0,0.4,0.3] );
hleg = legend('All','Axon guidance');  set(hleg,'Location','SouthEast'); legend('boxoff');     
set(gca,'box','off');  
set(gca,'ytick',[0:0.1:1]); set(gca,'xtick',[0:0.1:1]);
title('');



figure;
scatterScores(human6FlatResults, human_gene_info, brainspanFlatResults , developing_genes_info);
xlabel('human6 flat-tree scores', 'fontsize',20); ylabel('kang flat-tree scores', 'fontsize',20);

figure;
ploth = displayMultiDist({human6Results,human6RandomResults}, 'Triplets agreement score',{'All','Random'} ,[0.4,0.9]);
% ploth = displayMultiDist({human6Results,cahoy_human_subset.scores_neurons, cahoy_human_subset.scores_astro, cahoy_human_subset.scores_oligo,human6RandomResults}, 'Triplets agreement score',{'All','Neurons','Astrocytes','Oligodendrocytes','Random'} ,[0.4,0.9]);
% ploth = displayMultiDist({human6Results,housekeeping_subset.scores}, 'Triplets agreement score',{'All','House keeping'} ,[0.4,0.9]);
% ploth = displayMultiDist({human6Results,axon_guidance_subset.scores}, 'Triplets agreement score',{'All','Axon guidance'} ,[0.4,0.9]);
set(gca,'Fontsize',17);         %set(gca,'xscale','log');
set(gca,'ytick',[0:0.1:1]); set(gca,'xtick',[0:0.1:1]);
set(ploth,'LineWidth',  3 );

figure('name','human6');
displayDist(human6Results, human6RandomResults, 'Triplets agreement score'); %human6
set(gca,'ytick',[0:0.2:0.7]);
figure;
[qvaluesAboveRandHuman6, qvaluesBelowRandHuman6] = getEmpiricalPvalue(human6RandomResults, human6Results);

writeScoresInTable(human6Results, human_gene_info.gene_symbols, qvaluesAboveRandHuman6);

figure('name','kang');
displayDist(brainspanResults, brainspanRandomResults, 'Triplets agreement score'); %kang
set(gca,'ytick',[0:0.05:0.16]);
figure;
[qvaluesAboveRandBrainspan, qvaluesBelowRandBrainspan] = getEmpiricalPvalue(brainspanRandomResults, brainspanResults);

figure('name','flat score human6');
displayDist(human6FlatResults, human6FlatRandomResults, 'flat-tree score , human6');
figure;
[qvaluesAboveRandFlatHuman6, qvaluesBelowRandFlatHuman6] = getEmpiricalPvalue(human6FlatRandomResults, human6FlatResults);

figure('name','tree score kang');
displayDist(brainspanFlatResults, brainspanFlatRandomResults, 'tree-based score, kang');
figure;
[qvaluesAboveRandFlatBrainspan, qvaluesBelowRandFlatBrainspan] = getEmpiricalPvalue(brainspanFlatRandomResults, brainspanFlatResults);

figure('name','physical location, human6');
displayDist(human6PhysicalResults, human6PhysicalRandomResults, 'physical location, human6');
figure;
[qvaluesAboveRandPhysicalHuman6, qvaluesBelowRandPhysicalHuman6] = getEmpiricalPvalue(human6PhysicalRandomResults, human6PhysicalResults);

[scores, scores_rnd, gene_info] = compute_spearman('human6', []);
figure;
random_scores  = [scores_rnd{1};scores_rnd{2}; scores_rnd{3}];
displayDist(scores, random_scores, 'spearman-based score')


[brainspan_scores, brainspan_scores_rnd, ~] = compute_spearman('brainspan', []);
figure;
brainspan_scores_rnd  = [brainspan_scores_rnd{1};brainspan_scores_rnd{2}; brainspan_scores_rnd{3}];
displayDist(brainspan_scores, brainspan_scores_rnd, 'spearman-based score');

end

function writeScoresInTable(scores, geneNames, qvalues)
     [scores, sortInd] = sort(scores ,'descend');
     geneNames = geneNames(sortInd);
     qvalues = qvalues(sortInd);
     cellTogether = cat(2, cellstr( num2str(scores) ), geneNames, cellstr( num2str(qvalues) ) );
     writeCellToFile('treeGeneList.csv', cellTogether);

end
function scatterDots = scatterScores(scoresA, geneInfoA, scoresB, geneInfoB)
    [ind_for_A, ind_for_B] = mapGenes(geneInfoA.gene_symbols, geneInfoB.gene_symbols, geneInfoA.entrez_ids, geneInfoB.entrez_ids);
    intersect_scoresA = scoresA(ind_for_A);
    intersect_scoresB = scoresB(ind_for_B);
%     selectedGenesA = geneInfoA.gene_symbols(ind_for_A);
%     selectedGenesB = geneInfoB.gene_symbols(ind_for_B);
%     [sortedScores, sortInd] = sort(intersect_scoresA .* intersect_scoresB ,'descend');
%     sortedGenes = selectedGenesA(sortInd);

%     cellTogether = cat(2,num2cell(human6Results), human_gene_info.gene_symbols, num2cell(qvaluesAboveRandHuman6));
%     writeCellToFile('sortedTreeGenes.txt',sortedGenes);

    scatterDots = scatter(intersect_scoresA, intersect_scoresB,1,'filled');

    title(sprintf('correlation:  %g (pearson) %g (spearman)', corr(intersect_scoresA, intersect_scoresB), corr(intersect_scoresA, intersect_scoresB,'type','Spearman') ), 'fontsize',20);
%     set(gca,'yscale','log')
%     set(gca,'xscale','log')
end

function scatterSubsetScores(scoresA, geneInfoA, scoresB, geneInfoB, subset_entrez, subset_names,colorTriplet)
    [ind_for_A, ~] = mapGenes(geneInfoA.gene_symbols,subset_names, geneInfoA.entrez_ids, subset_entrez);
    [ind_for_B, ~] = mapGenes(geneInfoB.gene_symbols,subset_names, geneInfoB.entrez_ids, subset_entrez);
    
    scoresA = scoresA(ind_for_A);
    scoresB = scoresB(ind_for_B);
    geneInfoA.gene_symbols = geneInfoA.gene_symbols(ind_for_A);
    geneInfoA.entrez_ids = geneInfoA.entrez_ids(ind_for_A);
    geneInfoB.gene_symbols = geneInfoB.gene_symbols(ind_for_B);
    geneInfoB.entrez_ids = geneInfoB.entrez_ids(ind_for_B);
    
    scatterDots = scatterScores(scoresA, geneInfoA, scoresB, geneInfoB);
    
    set(scatterDots,'CData',colorTriplet);
    set(scatterDots,'SizeData',20);
    
end

function [scoresA, gene_symbols, entrez_ids] = addScoresToSubset(scoresA, geneInfo, subset_entrez, subset_names)
    [ind_for_A, ~] = mapGenes(geneInfo.gene_symbols,subset_names, geneInfo.entrez_ids, subset_entrez);
        
    scoresA = scoresA(ind_for_A);
    gene_symbols = geneInfo.gene_symbols(ind_for_A);
    entrez_ids = geneInfo.entrez_ids(ind_for_A);
end

function [qvaluesAboveRand, qvaluesBelowRand] = getEmpiricalPvalue(randomScores, scores)
    [numberOfSamples, repeats] = size(randomScores);
    assert(numberOfSamples == length(scores) );
    
    scoresRepeated = repmat(scores, [1, repeats]);
    pvaluesAboveRand = (1 + sum(scoresRepeated < randomScores,2))  ./ (repeats + 1);
    pvaluesBelowRand = (1 + sum(scoresRepeated > randomScores,2))  ./ (repeats + 1);
    
%     [~, qvaluesAboveRand] = mafdr(pvaluesAboveRand);
%     [~, qvaluesBelowRand] = mafdr(pvaluesBelowRand);
    [~, qvaluesAboveRand] = mafdr(pvaluesAboveRand,'lambda',0.15);
    [~, qvaluesBelowRand] = mafdr(pvaluesBelowRand,'lambda',0.15);
    
    hist(qvaluesAboveRand,100);
    title('significance to get a score above random','fontsize',20);
    xlabel('fdr corrected pvalues','fontsize',20);
    ylabel('number of genes','fontsize',20);
%     set(gca,'xscale','log')
end

function displayDist(scores, random_scores, xlabel_string)
    random_scores  = random_scores(:);
    min_val = min(scores);
    max_val = max(scores);

    spacing = linspace(min_val,max_val,100);
    vals_in_bins1 = histc(scores, spacing);
    vals_in_bins1 = vals_in_bins1 / sum(vals_in_bins1);
    vals_in_bins2 = histc(random_scores, spacing);
    vals_in_bins2 = vals_in_bins2 / sum(vals_in_bins2);

    together = [vals_in_bins1(:) , vals_in_bins2(:)];
    bar(spacing,together);
    xlabel(xlabel_string,'fontsize',20);
    h_legend = legend('Human genes','Random');
    legend('boxoff');
    set(h_legend,'FontSize',20);

    sortedRandom_scores = sort(random_scores);
    pointOf0001 = round( length(random_scores) * 0.001 ) ;
    top1000 = sortedRandom_scores(end-pointOf0001);
    bottom1000 = sortedRandom_scores(pointOf0001);

    countAbove1000 = sum( scores > top1000);
    countBelow1000 = sum( scores < bottom1000);

%     titleString = sprintf('%d genes < %g, %d genes > %g', countBelow1000, bottom1000, countAbove1000, top1000);
%     title(titleString,'fontsize',20);
    set(gca,'box','off');
    set(gca,'Fontsize',17);
end


function ploth = displayMultiDist(scores_cell, xlabel_string,legendStrings,range)
    spacing = linspace(range(1),range(2),25);
    vals_in_bins = nan(length(spacing), length(scores_cell));
    for i =1 :length(scores_cell);
        scores = scores_cell{i};
        vals_in_bins(:,i) = histc(scores(:), spacing);
        vals_in_bins(:,i) = vals_in_bins(:,i) / sum(vals_in_bins(:,i));
    end

    ploth = plot(spacing,vals_in_bins);
    xlabel(xlabel_string,'fontsize',20);
    h_legend = legend(legendStrings);
    legend('boxoff');
    set(h_legend,'FontSize',20);
    set(gca,'box','off');
    set(gca,'Fontsize',17);
    set(ploth,'LineWidth',  2 );
end