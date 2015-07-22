function viewScatterPlots()
init;
set(groot,'defaultAxesColorOrder','default');       

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
% colors = [all_samples_color; [ 0 0 1;  0 0 0.3; 1 0 0 ; 0.5,0.8,0.3] ];
% 
% set(0,'DefaultAxesColorOrder',colors);
% set(0,'DefaultAxesLineStyleOrder',{'-',':','--','.-'});

xtick = -0.1:0.1:1;
ytick = -0.1:0.1:1;

% xlimit = [0.4,0.9];
xlimit = [-0.13,0.8]; % full tree
% xlimit = [-0.1,0.3]; % for cortex
ylimit = [0, 0.25] ; % full tree
% ylimit = [0, 0.4] ; % for cortex


% load('data_matfile/treeResults_with_newKang100.mat');
%  load('data_matfile/simpleMeasurementTree.mat');
% load('data_matfile/fullTreeResultsTriplets.mat');
% load('data_matfile/simpleMeasurementFullTree.mat');
% load('data_matfile/cortexTreeResultsTriplets.mat');
% load('data_matfile/simpleMeasurementCortexTree.mat');
% load('data_matfile/cortexAllTreeResultsTriplets.mat');
% load('data_matfile/simpleMeasurementCortexAllTree.mat');
% load('data_matfile/simpleMeasurementAllRegionsCortexPar.mat');
% load('data_matfile/simpleMeasurementAllTreePar.mat');
% load('data_matfile/simpleMeasurementAllTreeNormPar.mat');

% load('data_matfile/simpleMeasurementAllRegionsPar.mat');
% load('data_matfile/simpleMeasurement16RegionsPar.mat');
% load('data_matfile/simpleMeasurementAllRegionsCortexPar.mat');
% load('data_matfile/simpleMeasurementCortexGrossPar.mat');

% data_results = load('results/human6GrossRegions-1-30000.mat');
data_results = load('results/human6AllRegions-1-30000.mat');
human6Results = data_results.results;
human6RandomResults = data_results.randomResults;
human_gene_info = data_results.gene_info;

% data_results = load('results/zapalaMouse-1-30000.mat');
data_results = load('results/kangAllRegions-1-30000.mat');
second_dataset_results = data_results.results;
second_dataset_random_results = data_results.randomResults;
second_dataset_genes_info = data_results.gene_info;

% homologous_genes_human6_zapala(human_gene_info, second_dataset_genes_info);

%load('data_matfile/kang_genes_2_fold.mat');
%error = drawVenn(getIndicesLargerThanThreshold(second_dataset_results, second_dataset_random_results, 0.01, 'right') , geneHasFoldChangeLargerThanX_dev_all_ages', correctedAnova_developing  < 0.01, 0.5);
%error = drawVenn(getIndicesLargerThanThreshold(human6Results, human6RandomResults, 0.01, 'right') , geneHasFoldChangeLargerThanX_human6', correctedAnova_human6  < 0.01, 0.5);
[~,sortInd] = sort(human6Results); medianIndex = sortInd(round(length(sortInd)/2) );
fprintf('%s ( %d ) is the median\n', human_gene_info.gene_symbols{medianIndex}, medianIndex);

%=============== counts how many genes are below 0.01 q-value============
% empirical_pvalue_human6 = getEmpiricalPvalues(human6Results, human6RandomResults);
% empirical_pvalue_human6 = mafdr(empirical_pvalue_human6, 'BHFDR', true);
% fprintf('all genes: %d / %d  (%d%%)\n', sum(empirical_pvalue_human6 < 0.01) , length(empirical_pvalue_human6) ,floor(sum(empirical_pvalue_human6< 0.01) / length(empirical_pvalue_human6)*100) );
% empirical_pvalue_developing = getEmpiricalPvalues(second_dataset_results, second_dataset_random_results);
% empirical_pvalue_developing = mafdr(empirical_pvalue_developing, 'BHFDR', true);
% fprintf('developing: %d / %d  (%d%%)\n', sum(empirical_pvalue_developing < 0.01) , length(empirical_pvalue_developing) ,floor(sum(empirical_pvalue_developing< 0.01) / length(empirical_pvalue_developing)*100) );
%========================================================================


figure('Name','qvalue vs gene count');
ploth = drawNumGenesVsRightPvalue(human6Results,human6RandomResults);
fileName = 'pvalue_cumsumm';
saveFigure(gcf, fileName, 'png');
saveFigure(gcf, fileName, 'tiff');
saveFigure(gcf, fileName, 'eps');
    
    

[cahoy_human_subset.scores_neurons, cahoy_human_subset.neurons_symbols, cahoy_human_subset.neuro_entrez] =...
    get_subset_scores('cahoy_neuro', human6Results, human_gene_info,human6RandomResults);
[cahoy_human_subset.scores_astro, cahoy_human_subset.astro_symbols, cahoy_human_subset.astro_entrez] =...
    get_subset_scores('cahoy_astro', human6Results, human_gene_info,human6RandomResults);
[cahoy_human_subset.scores_oligo, cahoy_human_subset.oligo_symbols, cahoy_human_subset.oligo_entrez] =...
    get_subset_scores('cahoy_oligo', human6Results, human_gene_info,human6RandomResults);
[housekeeping_subset.scores, housekeeping_subset.symbol, housekeeping_subset.entrez] = get_subset_scores(...
    'housekeeping', human6Results, human_gene_info, human6RandomResults);
[hoxGenes_subset.scores, hoxGenes_subset.symbol, hoxGenes_subset.entrez] = get_subset_scores(...
    'HOX', human6Results, human_gene_info, human6RandomResults);
[axon_guidance_subset.scores, axon_guidance_subset.symbol, axon_guidance_subset.entrez] = get_subset_scores(...
    'axon_guidance', human6Results, human_gene_info,human6RandomResults);
[paxGene.scores, paxGene.symbol, paxGene.entrez] = get_subset_scores(...
    'PAX', human6Results, human_gene_info,human6RandomResults);
[soxGene.scores, soxGene.symbol, soxGene.entrez] = get_subset_scores(...
    'SOX2', human6Results, human_gene_info,human6RandomResults);
[serotoninGene.scores, serotoninGene.symbols, serotoninGene.entrez] = get_subset_scores(...
    'serotonin', human6Results, human_gene_info,human6RandomResults);
[dopaminGene.scores, dopaminGene.symbols, dopaminGene.entrez] = get_subset_scores(...
    'dopamin', human6Results, human_gene_info,human6RandomResults);
[dopaminAndSertonin.scores, dopaminAndSertonin.symbols, dopaminAndSertonin.entrez] = get_subset_scores(...
    'dopaminAndSertonin', human6Results, human_gene_info,human6RandomResults);
age_scores = {}; age_labels = {};
for i =5
    [age_scores_i,~,~,age_labels_i] = get_subset_scores(...
    sprintf('age-%d',i ), human6Results, human_gene_info,human6RandomResults);
    age_scores = cat(2,age_scores, {age_scores_i});
    age_labels = cat(2,age_labels, {age_labels_i});
end
age_colors = autumn(length(age_labels) +1);
age_colors = age_colors(1:end-1,:);


fprintf('human6: ');
moreThenXPercent(human6Results, human6RandomResults, 0.99);

fprintf('second dataset: ');
moreThenXPercent(second_dataset_results, second_dataset_random_results, 0.99);



% figure('name','Scatter with random');
% hold on;
% scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
% randomScatterDots = scatterScores(human6RandomResults, human_gene_info, second_dataset_random_results , second_dataset_genes_info);
% set(randomScatterDots, 'CData',random_samples_color);
% xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
% hleg = legend('All','Random');  set(hleg,'Location','Northwest'); legend('boxoff');    set(hleg,'FontSize',20); 
% set(gca,'box','off');  
% set(gca,'ytick',ytick); set(gca,'xtick',xtick);
% title('');
% redoScatterImage(xlimit);
% saveFigure(gcf, 'figures/scatter_random.png', 'png');
% saveFigure(gcf, 'figures/scatter_random.tiff', 'tiff');
% saveFigure(gcf, 'figures/scatter_random', 'eps');
% 
% figure('name','Cell type scatter');
% hold on;
% scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, cahoy_human.oligo_human_entrez, cahoy_human.oligo_human_hsbc, oligo_color);
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, cahoy_human.astro_human_entrez, cahoy_human.astro_human_hsbc, astro_color);
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, cahoy_human.neurons_human_entrez, cahoy_human.neurons_human_hsbc, neurons_color);
% 
% xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
% hleg = legend('All','Oligodendrocytes','Astrocytes','Neurons');  set(hleg,'Location','Northwest'); legend('boxoff');    set(hleg,'FontSize',20); 
% set(gca,'box','off');  
% set(gca,'ytick',ytick); set(gca,'xtick',xtick);
% title('');
% redoScatterImage(xlimit);
% saveFigure(gcf, 'figures/scatter_cellType.png', 'png');
% saveFigure(gcf, 'figures/scatter_cellType.tiff', 'tiff');
% saveFigure(gcf, 'figures/scatter_cellType.eps', 'eps');
% 
% figure('name','Housekeeping scatter');
% scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
% xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
% hold on;
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, housekeeping.entrez ,housekeeping.symbol, house_keeping_color );
% hleg = legend('All','Housekeeping');  set(hleg,'Location','Northwest'); legend('boxoff');  set(hleg,'FontSize',20);   
% set(gca,'box','off');  
% set(gca,'ytick',ytick); set(gca,'xtick',xtick);
% title('');
% redoScatterImage(xlimit);
% saveFigure(gcf, 'figures/scatter_house.png', 'png');
% saveFigure(gcf, 'figures/scatter_house.tiff', 'tiff');
% saveFigure(gcf, 'figures/scatter_house.eps', 'eps');
% 
% figure('name','Hox and Axon guidance scatter');
% scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
% set(scatterDots,'SizeData',40);
% xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
% hold on;
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan(size(axon_guidance.symbol )) ,axon_guidance.symbol , axon_guidance_color );
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan(size(hoxGenes.symbol )) ,hoxGenes.symbol , hox_color );
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan ,paxGeneSymbol , pax_color );
% % scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, soxGene.entrez ,{'SOX2'} , [1 0 0] );
% hleg = legend('All','Axon guidance','Hox family','Pax family');  set(hleg,'Location','Northwest'); legend('boxoff');     set(hleg,'FontSize',20);
% set(gca,'box','off');  
% set(gca,'ytick',ytick); set(gca,'xtick',xtick);
% title('');
% redoScatterImage(xlimit);
% saveFigure(gcf, 'figures/scatter_axon.png', 'png');
% saveFigure(gcf, 'figures/scatter_axon.tiff', 'tiff');
% saveFigure(gcf, 'figures/scatter_axon', 'eps');
% 
% 
% 
% figure('name','Serotonin scatter');
% scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
% set(scatterDots,'SizeData',40);
% xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement (Kang-2011)', 'fontsize',20);
% hold on;
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan ,dopaminAndSertoninGeneSymbol , pax_color ,false);
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan(size(serotoninGeneSymbol )) ,serotoninGeneSymbol , axon_guidance_color ,false);
% scatterSubsetScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info, nan(size(dopaminGeneSymbol )) ,dopaminGeneSymbol , hox_color ,false);
% hleg = legend('All','Serotonin & Dopamine', 'Serotonin','Dopamine');  set(hleg,'Location','Northwest'); legend('boxoff');     set(hleg,'FontSize',20);
% set(gca,'box','off');  
% set(gca,'ytick',ytick); set(gca,'xtick',xtick);
% title('');
% redoScatterImage(xlimit);
% saveFigure(gcf, 'figures/scatter_serotonin.png', 'png');
% saveFigure(gcf, 'figures/scatter_serotonin.tiff', 'tiff');
% saveFigure(gcf, 'figures/scatter_serotonin', 'eps');
% 


% 
% 
% figure;
% scatterScores(human6FlatResults, human_gene_info, brainspanFlatResults , second_dataset_genes_info);
% xlabel('human6 flat-tree scores', 'fontsize',20); ylabel('kang flat-tree scores', 'fontsize',20);


% compareDistributions(human6Results, human6RandomResults(:), 'all human6', 'randomHuman6');
% compareDistributions( second_dataset_results,second_dataset_random_results(:), 'all kang','randomKang');


draw_groups_dist('Human6',{human6Results,human6RandomResults}, ...
    {'All','Random'}, [all_samples_color;random_samples_color], xlimit,ylimit);
draw_groups_dist('Kang',{second_dataset_results,second_dataset_random_results}, ...
    {'All','Random'}, [all_samples_color;random_samples_color], xlimit,ylimit);
draw_groups_dist('Human6 cahoy',{  human6Results, cahoy_human_subset.scores_oligo, cahoy_human_subset.scores_astro, cahoy_human_subset.scores_neurons},...
    {'All','Oligodendrocytes','Astrocytes','Neurons'} , [all_samples_color;oligo_color;astro_color;neurons_color], xlimit,ylimit);
draw_groups_dist('Human6 housekeeping',{human6Results,housekeeping_subset.scores}, ...
    {'All','Housekeeping'} , [all_samples_color;house_keeping_color], xlimit,ylimit);
draw_groups_dist('Human6 axon guidance',{human6Results,axon_guidance_subset.scores,hoxGenes_subset.scores, paxGene.scores},...
    {'All','Axon guidance','Hox','Pax'}, [all_samples_color;axon_guidance_color;hox_color;pax_color], xlimit,ylimit);
draw_groups_dist('Human6 serotonin',{human6Results,serotoninGene.scores}, ...
    {'All','Serotonin'}, [all_samples_color;pax_color], xlimit,ylimit);
draw_groups_dist('Age',[{human6Results}, age_scores], ...
    [{'All'}, age_labels], [all_samples_color;age_colors], xlimit,ylimit);

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

