function viewScatterPlots()

random_samples_color = [0.4, 0.4,0.4];
all_samples_color = [0.8, 0.8,0.8];
house_keeping_color = [0.4,0.7,0] ;
hox_color = [.7 .9 .1];
axon_guidance_color =  [ 0.7, 0.4, 0.9];
pax_color = [ 0.9, 0.2, 0.5];

neurons_color =  [ 0, 0, 1];
oligo_color =  [ 0,0,0.3];
astro_color =  [ 1, 0, 0];

% colors = distinguishable_colors(20);
colors = [all_samples_color; [ 0 0 1;  0 0 0.3; 1 0 0 ; 0.5,0.8,0.3] ];

set(0,'DefaultAxesColorOrder',colors);
set(0,'DefaultAxesLineStyleOrder',{'-',':','--','.-'});

xtick = -0.1:0.1:1;
ytick = -0.1:0.1:1;

% xlimit = [0.4,0.9];
xlimit = [-0.13,0.8]; % full tree
% xlimit = [-0.1,0.3]; % for cortex
ylimit = [0, 0.25] ; % full tree
% ylimit = [0, 0.4] ; % for cortex


% load('treeResults_with_newKang100.mat');
%  load('simpleMeasurementTree.mat');
% load('fullTreeResultsTriplets.mat');
% load('simpleMeasurementFullTree.mat');
% load('cortexTreeResultsTriplets.mat');
% load('simpleMeasurementCortexTree.mat');
% load('cortexAllTreeResultsTriplets.mat');
% load('simpleMeasurementCortexAllTree.mat');
% load('simpleMeasurementAllRegionsCortexPar.mat');
% load('simpleMeasurementAllTreePar.mat');
% load('simpleMeasurementAllTreeNormPar.mat');

load('simpleMeasurementAllRegionsPar.mat');
% load('simpleMeasurement16RegionsPar.mat');
% load('simpleMeasurementAllRegionsCortexPar.mat');
% load('simpleMeasurementCortexGrossPar.mat');


%load('kang_genes_2_fold.mat');
%error = drawVenn(getIndicesLargerThanThreshold(brainspanResults, brainspanRandomResults, 0.01, 'right') , geneHasFoldChangeLargerThanX_dev_all_ages', correctedAnova_developing  < 0.01, 0.5);
%error = drawVenn(getIndicesLargerThanThreshold(human6Results, human6RandomResults, 0.01, 'right') , geneHasFoldChangeLargerThanX_human6', correctedAnova_human6  < 0.01, 0.5);
[~,sortInd] = sort(human6Results); medianIndex = sortInd(round(length(sortInd)/2) );
fprintf('%s ( %d ) is the median\n', human_gene_info.gene_symbols{medianIndex}, medianIndex);

%=============== counts how many genes are below 0.01 q-value============
% empirical_pvalue_human6 = getEmpiricalPvalues(human6Results, human6RandomResults);
% empirical_pvalue_human6 = mafdr(empirical_pvalue_human6, 'BHFDR', true);
% fprintf('all genes: %d / %d  (%d%%)\n', sum(empirical_pvalue_human6 < 0.01) , length(empirical_pvalue_human6) ,floor(sum(empirical_pvalue_human6< 0.01) / length(empirical_pvalue_human6)*100) );
% empirical_pvalue_developing = getEmpiricalPvalues(brainspanResults, brainspanRandomResults);
% empirical_pvalue_developing = mafdr(empirical_pvalue_developing, 'BHFDR', true);
% fprintf('developing: %d / %d  (%d%%)\n', sum(empirical_pvalue_developing < 0.01) , length(empirical_pvalue_developing) ,floor(sum(empirical_pvalue_developing< 0.01) / length(empirical_pvalue_developing)*100) );
%========================================================================


figure('Name','qvalue vs gene count');
ploth = drawNumGenesVsRightPvalue(human6Results,human6RandomResults);
%     set(gca,'ytick',[0:5000:20000]);
fileName = 'pvalue_cumsumm';
saveFigure(gcf, fileName, 'png');
saveFigure(gcf, fileName, 'tiff');
saveFigure(gcf, fileName, 'eps');
    
    
    
cahoy_human = load('cahoy_human_homolgues.mat');
[cahoy_human_subset.scores_neurons, cahoy_human_subset.neurons_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.neurons_human_entrez, cahoy_human.neurons_human_hsbc);
[cahoy_human_subset.scores_astro, cahoy_human_subset.astro_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.astro_human_entrez, cahoy_human.astro_human_hsbc);
[cahoy_human_subset.scores_oligo, cahoy_human_subset.oligo_symbols, cahoy_human_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, cahoy_human.oligo_human_entrez, cahoy_human.oligo_human_hsbc);
cahoy_human_subset.empirical_pvalue_neurons = getEmpiricalPvalues(cahoy_human_subset.scores_neurons, human6RandomResults);
cahoy_human_subset.empirical_pvalue_neurons = mafdr(cahoy_human_subset.empirical_pvalue_neurons, 'BHFDR', true);
cahoy_human_subset.empirical_pvalue_astro = getEmpiricalPvalues(cahoy_human_subset.scores_astro, human6RandomResults);
cahoy_human_subset.empirical_pvalue_astro = mafdr(cahoy_human_subset.empirical_pvalue_astro, 'BHFDR', true);
cahoy_human_subset.empirical_pvalue_oligo = getEmpiricalPvalues(cahoy_human_subset.scores_oligo, human6RandomResults);
cahoy_human_subset.empirical_pvalue_oligo = mafdr(cahoy_human_subset.empirical_pvalue_oligo, 'BHFDR', true);
fprintf('neurons: %d / %d  (%d%%)\n', sum(cahoy_human_subset.empirical_pvalue_neurons < 0.01) , length(cahoy_human_subset.empirical_pvalue_neurons) ,floor(sum(cahoy_human_subset.empirical_pvalue_neurons < 0.01) / length(cahoy_human_subset.empirical_pvalue_neurons)*100) );
fprintf('astro: %d / %d  (%d%%)\n', sum(cahoy_human_subset.empirical_pvalue_astro < 0.01) , length(cahoy_human_subset.empirical_pvalue_astro) ,floor(sum(cahoy_human_subset.empirical_pvalue_astro < 0.01) / length(cahoy_human_subset.empirical_pvalue_astro)*100) );
fprintf('oligo: %d / %d  (%d%%)\n', sum(cahoy_human_subset.empirical_pvalue_oligo < 0.01) , length(cahoy_human_subset.empirical_pvalue_oligo) ,floor(sum(cahoy_human_subset.empirical_pvalue_oligo < 0.01) / length(cahoy_human_subset.empirical_pvalue_oligo)*100) );



[housekeeping.ensembl_gene_id,housekeeping.transcript_gene_id ,housekeeping.entrez ,housekeeping.symbol ] = textread('/cortex/data/gene_sets/house_keeping_genes/Erez_2003.txt','%s %s %d %s','headerlines',1);
[housekeeping_subset.scores, housekeeping_subset.symbols, housekeeping_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, housekeeping.entrez, housekeeping.symbol);

hoxGenesSymbols = {'Hoxa1';'Hoxa10';'Hoxa11';'Hoxa2';'Hoxa3';'Hoxa4';'Hoxa5';'Hoxa6';'Hoxa7';'Hoxa9';'Hoxb1';'Hoxb13';'Hoxb3';'Hoxb4';'Hoxb5';'Hoxb6';'Hoxb9';'Hoxc10';'Hoxc12';'Hoxc13';'Hoxc4';'Hoxc5';'Hoxc8';'Hoxc9';'Hoxd1';'Hoxd12';'Hoxd13';'Hoxd3';'Hoxd4';'Hoxd8';'Hoxd8';'Hoxd9';};
[~, indInList] = ismember(upper(hoxGenesSymbols),human_gene_info.gene_symbols);
hoxGenes.symbol = human_gene_info.gene_symbols(indInList);
hoxGenes.entrez = human_gene_info.entrez_ids(indInList);
[hoxGenes_subset.scores, hoxGenes_subset.symbols, hoxGenes_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, hoxGenes.entrez, hoxGenes.symbol);
hoxGenes_subset.empirical_pvalue = getEmpiricalPvalues(hoxGenes_subset.scores, human6RandomResults);
hoxGenes_subset.empirical_pvalue = mafdr(hoxGenes_subset.empirical_pvalue, 'BHFDR', true);
fprintf('hox: %d / %d  (%d%%)\n', sum(hoxGenes_subset.empirical_pvalue < 0.01) , length(hoxGenes_subset.empirical_pvalue) ,floor(sum(hoxGenes_subset.empirical_pvalue < 0.01) / length(hoxGenes_subset.empirical_pvalue)*100) );

[axon_guidance.symbol ] = textread('/cortex/data/gene_sets/human_axon_gudiance','%s');
[axon_guidance_subset.scores, axon_guidance_subset.symbols, axon_guidance_subset.entrez] = addScoresToSubset(human6Results, human_gene_info, nan(size(axon_guidance.symbol)), axon_guidance.symbol);


paxGeneSymbol = {'Pax1','Pax2','Pax3','Pax4','Pax5','Pax6','Pax7','Pax8';};
[~, indInList] = ismember(upper(paxGeneSymbol),human_gene_info.gene_symbols);
paxGene.symbol = human_gene_info.gene_symbols(indInList);
paxGene.entrez = human_gene_info.entrez_ids(indInList);
[paxGene.scores, paxGene.symbols, paxGene.entrez] = addScoresToSubset(human6Results, human_gene_info, paxGene.entrez, paxGene.symbol);

sox2GeneSymbol = {'SOX2';};
[~, indInList] = ismember(upper(sox2GeneSymbol),human_gene_info.gene_symbols);
soxGene.symbol = human_gene_info.gene_symbols(indInList);
soxGene.entrez = human_gene_info.entrez_ids(indInList);
[soxGene.scores, soxGene.symbols, soxGene.entrez] = addScoresToSubset(human6Results, human_gene_info, soxGene.entrez, soxGene.symbol);

moreThenXPercent(human6Results, human6RandomResults, 0.99);



serotoninGeneSymbol = {'HTR1A','HTR1B','HTR1D','HTR1E','HTR1F','HTR2A','HTR2B','HTR2C','HTR3A','HTR3B','HTR3C','HTR3D','HTR3E','HTR4','HTR5A','HTR6','HTR7','TPH1','TPH2', 'SLC6A4'};
[~, indInList] = ismember(upper(serotoninGeneSymbol),human_gene_info.gene_symbols);
serotoninGene.symbol = human_gene_info.gene_symbols(indInList);
serotoninGene.entrez = human_gene_info.entrez_ids(indInList);
[serotoninGene.scores, serotoninGene.symbols, serotoninGene.entrez] = addScoresToSubset(human6Results, human_gene_info, serotoninGene.entrez, serotoninGene.symbol);

dopaminGeneSymbol = {'COMT','DRD1','DRD2','DRD3','DRD4','DRD5','SLC6A3','TH'};
[~, indInList] = ismember(upper(dopaminGeneSymbol),human_gene_info.gene_symbols);
dopaminGene.symbol = human_gene_info.gene_symbols(indInList);
dopaminGene.entrez = human_gene_info.entrez_ids(indInList);
[dopaminGene.scores, dopaminGene.symbols, dopaminGene.entrez] = addScoresToSubset(human6Results, human_gene_info, dopaminGene.entrez, dopaminGene.symbol);

dopaminAndSertoninGeneSymbol = {'DDC', 'MAOA','MAOB','SLC18A1','SLC18A2'};
[~, indInList] = ismember(upper(dopaminAndSertoninGeneSymbol),human_gene_info.gene_symbols);
dopaminAndSertoninGene.symbol = human_gene_info.gene_symbols(indInList);
dopaminAndSertoninGene.entrez = human_gene_info.entrez_ids(indInList);
[dopaminAndSertoninGene.scores, dopaminAndSertoninGene.symbols, dopaminAndSertoninGene.entrez] = addScoresToSubset(human6Results, human_gene_info, dopaminAndSertoninGene.entrez, dopaminAndSertoninGene.symbol);




figure('name','Scatter with random');
hold on;
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
randomScatterDots = scatterScores(human6RandomResults, human_gene_info, brainspanRandomResults , developing_genes_info);
set(randomScatterDots, 'CData',random_samples_color);
xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
hleg = legend('All','Random');  set(hleg,'Location','Northwest'); legend('boxoff');    set(hleg,'FontSize',20); 
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'scatter_random.png', 'png');
saveFigure(gcf, 'scatter_random.tiff', 'tiff');
saveFigure(gcf, 'scatter_random', 'eps');

figure('name','Cell type scatter');
hold on;
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.oligo_human_entrez, cahoy_human.oligo_human_hsbc, oligo_color);
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.astro_human_entrez, cahoy_human.astro_human_hsbc, astro_color);
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, cahoy_human.neurons_human_entrez, cahoy_human.neurons_human_hsbc, neurons_color);

xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
hleg = legend('All','Oligodendrocytes','Astrocytes','Neurons');  set(hleg,'Location','Northwest'); legend('boxoff');    set(hleg,'FontSize',20); 
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'scatter_cellType.png', 'png');
saveFigure(gcf, 'scatter_cellType.tiff', 'tiff');
saveFigure(gcf, 'scatter_cellType.eps', 'eps');

figure('name','Housekeeping scatter');
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, housekeeping.entrez ,housekeeping.symbol, house_keeping_color );
hleg = legend('All','Housekeeping');  set(hleg,'Location','Northwest'); legend('boxoff');  set(hleg,'FontSize',20);   
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'scatter_house.png', 'png');
saveFigure(gcf, 'scatter_house.tiff', 'tiff');
saveFigure(gcf, 'scatter_house.eps', 'eps');

figure('name','Hox and Axon guidance scatter');
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
set(scatterDots,'SizeData',40);
xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan(size(axon_guidance.symbol )) ,axon_guidance.symbol , axon_guidance_color );
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan(size(hoxGenes.symbol )) ,hoxGenes.symbol , hox_color );
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan ,paxGeneSymbol , pax_color );
% scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, soxGene.entrez ,{'SOX2'} , [1 0 0] );
hleg = legend('All','Axon guidance','Hox family','Pax family');  set(hleg,'Location','Northwest'); legend('boxoff');     set(hleg,'FontSize',20);
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'scatter_axon.png', 'png');
saveFigure(gcf, 'scatter_axon.tiff', 'tiff');
saveFigure(gcf, 'scatter_axon', 'eps');



figure('name','Serotonin scatter');
scatterDots = scatterScores(human6Results, human_gene_info, brainspanResults , developing_genes_info);
set(scatterDots,'SizeData',40);
xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
hold on;
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan ,dopaminAndSertoninGeneSymbol , pax_color ,false);
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan(size(serotoninGeneSymbol )) ,serotoninGeneSymbol , axon_guidance_color ,false);
scatterSubsetScores(human6Results, human_gene_info, brainspanResults , developing_genes_info, nan(size(dopaminGeneSymbol )) ,dopaminGeneSymbol , hox_color ,false);
hleg = legend('All','Serotonin & Dopamine', 'Serotonin','Dopamine');  set(hleg,'Location','Northwest'); legend('boxoff');     set(hleg,'FontSize',20);
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'scatter_serotonin.png', 'png');
saveFigure(gcf, 'scatter_serotonin.tiff', 'tiff');
saveFigure(gcf, 'scatter_serotonin', 'eps');



% 
% 
% figure;
% scatterScores(human6FlatResults, human_gene_info, brainspanFlatResults , developing_genes_info);
% xlabel('human6 flat-tree scores', 'fontsize',20); ylabel('kang flat-tree scores', 'fontsize',20);


% compareDistributions(human6Results, human6RandomResults(:), 'all human6', 'randomHuman6');
% compareDistributions( brainspanResults,brainspanRandomResults(:), 'all kang','randomKang');
compareDistributions(cahoy_human_subset.scores_neurons, human6Results, 'neurons','all human6');
compareDistributions(cahoy_human_subset.scores_oligo, human6Results,   'oligo', 'all human6');
compareDistributions(cahoy_human_subset.scores_astro, human6Results,  'astro', 'all human6');
compareDistributions(axon_guidance_subset.scores,human6Results,  'axon', 'all human6');
compareDistributions( hoxGenes_subset.scores, human6Results, 'hox','all human6');
compareDistributions( paxGene.scores, human6Results, 'pax','all human6');
compareDistributions(housekeeping_subset.scores, human6Results,  'house', 'all human6');

% compareDistributions( hoxGenes_subset.scores, human6RandomResults(:), 'hox','randomHuman6');
% compareDistributions( housekeeping_subset.scores, human6RandomResults(:), 'house','randomHuman6');

f = figure('name','Human6 distribution');  
ploth = displayMultiDist({human6Results,human6RandomResults}, 'BRO-agreement (ABA6-2013)',{'All','Random'} ,xlimit);
set(ploth(2), 'Color', random_samples_color );  set(ploth,'LineWidth',  3 );
redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
set(gca,'ytick',[0:0.2:1]);
saveFigure(gcf, 'dist_random.png', 'png');
saveFigure(gcf, 'dist_random.tiff', 'tiff');
saveFigure(gcf, 'dist_random', 'eps');

figure('name','Kang distribution');  
ploth = displayMultiDist({brainspanResults,brainspanRandomResults}, 'BRO-agreement (Kang-2011)',{'All','Random'} ,xlimit);
set(ploth(2), 'Color', random_samples_color );  set(ploth,'LineWidth',  3 );
redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
set(gca,'ytick',[0:0.2:1]);
saveFigure(gcf, 'dist_kang_random.png', 'png');
saveFigure(gcf, 'dist_kang_random.tiff', 'tiff');
saveFigure(gcf, 'dist_kang_random', 'eps');

figure('name','Human6 cahoy distribution'); 
ploth = displayMultiDist({  human6Results, cahoy_human_subset.scores_oligo, cahoy_human_subset.scores_astro, cahoy_human_subset.scores_neurons}, 'BRO-agreement (ABA6-2013)',{'All','Oligodendrocytes','Astrocytes','Neurons'} ,xlimit);
set(ploth(2),'Color', oligo_color ); 
set(ploth(3),'Color', astro_color ); 
set(ploth(4),'Color', neurons_color ); 
redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
ylim(ylimit);
saveFigure(gcf, 'dist_cellType.png', 'png');
saveFigure(gcf, 'dist_cellType.tiff', 'tiff');
saveFigure(gcf, 'dist_cellType', 'eps');

figure('name','Human6 housekeeping distribution'); 
ploth = displayMultiDist({human6Results,housekeeping_subset.scores}, 'BRO-agreement (ABA6-2013)',{'All','Housekeeping'} ,xlimit);
set(ploth(2),'Color', house_keeping_color );  set(ploth,'LineWidth',  3 );
redoDistImage(xlimit);
ylim(ylimit);
saveFigure(gcf, 'dist_house.png', 'png');
saveFigure(gcf, 'dist_house.tiff', 'tiff');
saveFigure(gcf, 'dist_house', 'eps');

figure('name','Human6 axon guidance distribution'); 
ploth = displayMultiDist({human6Results,axon_guidance_subset.scores,hoxGenes_subset.scores, paxGene.scores}, 'BRO-agreement (ABA6-2013)',{'All','Axon guidance','Hox','Pax'} ,xlimit);
set(ploth(2),'Color', axon_guidance_color ); 
set(ploth(3),'Color', hox_color ); 
set(ploth(4),'Color', pax_color ); 
set(ploth,'LineWidth',  3 );
redoDistImage(xlimit);
ylim(ylimit);
saveFigure(gcf, 'dist_axon.png', 'png');
saveFigure(gcf, 'dist_axon.tiff', 'tiff');
saveFigure(gcf, 'dist_axon', 'eps');


figure('name','Human6 serotonin distribution'); 
ploth = displayMultiDist({human6Results,serotoninGene.scores}, 'BRO-agreement (ABA6-2013)',{'All','Serotonin'} ,xlimit);
set(ploth(2), 'Color', pax_color ); 
set(ploth,'LineWidth',  3 );
redoDistImage(xlimit);
ylim(ylimit);
saveFigure(gcf, 'dist_serotonin.png', 'png');
saveFigure(gcf, 'dist_serotonin.tiff', 'tiff');
saveFigure(gcf, 'dist_serotonin', 'eps');


end

function significant_scores = getIndicesLargerThanThreshold(scores, randomScores, p_value, tail)
    randomScores = randomScores(:);
    scores = scores(:);
    
    switch tail
        case 'right'
            sortedRandomScores = sort(randomScores);
            indexAtThreshold = length(randomScores) - round(p_value * length(randomScores));
            valueAtThreshold = sortedRandomScores(indexAtThreshold);
            significant_scores = scores > valueAtThreshold;
        case 'left'
            sortedRandomScores = sort(randomScores,'descend');
            indexAtThreshold = length(randomScores) - round(p_value * length(randomScores));
            valueAtThreshold = sortedRandomScores(indexAtThreshold);
            significant_scores = scores < valueAtThreshold;
        case 'both'
            significant_scores_left = getIndicesLargerThanThreshold(scores, randomScores, 0.5*p_value, 'left');
            significant_scores_right = getIndicesLargerThanThreshold(scores, randomScores, 0.5*p_value, 'right');
            significant_scores = significant_scores_left | significant_scores_right;
    end
            
end

function scatterDots =  scatterScores(scoresA, geneInfoA, scoresB, geneInfoB, showNames)
    if ~exist('showNames','var')
        showNames = false;
    end
    [ind_for_A, ind_for_B] = mapGenes(geneInfoA.gene_symbols, geneInfoB.gene_symbols, geneInfoA.entrez_ids, geneInfoB.entrez_ids);
    intersect_scoresA = scoresA(ind_for_A);
    intersect_scoresB = scoresB(ind_for_B);
    scatterDots = scatter(intersect_scoresA, intersect_scoresB,1,'filled');
    if (showNames)
        intersect_symbols = geneInfoA.gene_symbols(ind_for_A);
        text( double(intersect_scoresA) ,double(intersect_scoresB) , intersect_symbols, 'horizontal','left', 'vertical','bottom','fontsize',10);
    end
%     scatterDots = transperntScatter(intersect_scoresA,intersect_scoresB,0.01,0.1);
    title(sprintf('correlation:  %g (pearson) %g (spearman)', corr(intersect_scoresA, intersect_scoresB), corr(intersect_scoresA, intersect_scoresB,'type','Spearman') ), 'fontsize',20);
    set(scatterDots,'SizeData',50);
end

function compareDistributions(distA, distB, nameA, nameB)
    pvalue = ranksum(distA, distB);
    meadianA = median(distA);
    meadianB = median(distB);
    fprintf('%s (median %g) - %s (median %g):\t\t %g both', nameA,meadianA,nameB,meadianB, pvalue);
    pvalueright = ranksum(distA, distB,'tail','right');
    fprintf(',\t\t %g right', pvalueright);
    pvalueleft = ranksum(distA, distB,'tail','left');
    fprintf(',\t\t %g left', pvalueleft);
    fprintf('\n');
 
end

function ploth = drawNumGenesVsRightPvalue(scores, randomScores)

 
    empirical_pvalue_human6 = getEmpiricalPvalues(scores, randomScores(:) );
    empirical_qvalue_human6 = mafdr(empirical_pvalue_human6, 'BHFDR', true);
   
    pvalueright = sort(empirical_qvalue_human6,'descend');
    ploth = plot(pvalueright, length(pvalueright):-1:1);
    set(ploth,'LineWidth',  3 );
    xlabel('q-value (FDR)');
    ylabel('number of genes');
    set(gca,'Fontsize',20);         
    set(gca,'xscale','log');
   
%     set(gca,'xscale','log');
%     xtick = -0.2:0.2:1;
%     ytick = -0.2:0.1:1;
%     set(gca,'ytick',ytick); set(gca,'xtick',xtick);
end

function showMedianDist(scores, expressionDistances,treeDistances)
    [~, medianInd] = sort(scores);
    medianGeneExpression = 3;
end

function moreThenXPercent(scores, randomScores, x_percent)

randomScoreSorted = sort(randomScores(:) );

index = round( x_percent * length(randomScoreSorted) );
how_many_larger = sum(scores > randomScoreSorted(index) );

fprintf( '%d%% of the genes are larger (%d)\n', round( how_many_larger / length(scores)*100), how_many_larger);
end

function redoDistImage(xlimit)
    xtick = -0.2:0.2:1;
    ytick = -0.2:0.1:1;

    set(gca,'Fontsize',20);         %set(gca,'xscale','log');
    set(gca,'ytick',ytick); set(gca,'xtick',xtick);
    xlim( xlimit)
end

function redoScatterImage(limit)
    xtick = -0.2:0.2:1;
    ytick = -0.2:0.2:1;

    set(gca,'Fontsize',20);         %set(gca,'xscale','log');
    set(gca,'ytick',ytick); set(gca,'xtick',xtick);
    xlim( limit )
    ylim( limit )
%     axis equal
end

function scatterDots = scatterSubsetScores(scoresA, geneInfoA, scoresB, geneInfoB, subset_entrez, subset_names,colorTriplet, showNames)
    if ~exist('showNames','var')
        showNames = false;
    end
    [ind_for_A, ~] = mapGenes(geneInfoA.gene_symbols,subset_names, geneInfoA.entrez_ids, subset_entrez);
    [ind_for_B, ~] = mapGenes(geneInfoB.gene_symbols,subset_names, geneInfoB.entrez_ids, subset_entrez);
    
    scoresA = scoresA(ind_for_A);
    scoresB = scoresB(ind_for_B);
    geneInfoA.gene_symbols = geneInfoA.gene_symbols(ind_for_A);
    geneInfoA.entrez_ids = geneInfoA.entrez_ids(ind_for_A);
    geneInfoB.gene_symbols = geneInfoB.gene_symbols(ind_for_B);
    geneInfoB.entrez_ids = geneInfoB.entrez_ids(ind_for_B);
    
    scatterDots = scatterScores(scoresA, geneInfoA, scoresB, geneInfoB, showNames);
    set(scatterDots,'CData',colorTriplet);
%      set(scatterDots,'FaceColor',colorTriplet);
    
    
end

function [scoresA, gene_symbols, entrez_ids] = addScoresToSubset(scoresA, geneInfo, subset_entrez, subset_names)
    [ind_for_A, ~] = mapGenes(geneInfo.gene_symbols,subset_names, geneInfo.entrez_ids, subset_entrez);
        
    scoresA = scoresA(ind_for_A);
    gene_symbols = geneInfo.gene_symbols(ind_for_A);
    entrez_ids = geneInfo.entrez_ids(ind_for_A);
end

function ploth = displayMultiDist(scores_cell, xlabel_string,legendStrings,range)
    spacing = linspace(range(1),range(2),50);
    vals_in_bins = nan(length(spacing), length(scores_cell));
    for i =1 :length(scores_cell);
        scores = scores_cell{i};
        vals_in_bins(:,i) = histc(scores(:), spacing);
%         vals_in_bins(:,i) = smooth(vals_in_bins(:,i), 0.1, 'moving');
        vals_in_bins(:,i) = vals_in_bins(:,i) / sum(vals_in_bins(:,i));
    end

    vals_in_bins_smooth = smooth(vals_in_bins,0.02);
    vals_in_bins_smooth = reshape(vals_in_bins_smooth,size(vals_in_bins));
    
    ploth = plot(spacing,  vals_in_bins_smooth    );
    xlabel(xlabel_string,'fontsize',20);
    h_legend = legend(legendStrings);
    legend('boxoff');
    set(h_legend,'FontSize',20);
    set(gca,'box','off');
    set(gca,'Fontsize',17);
    set(ploth,'LineWidth',  2 );
end