function viewScatterPlots_human_mouse()

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


% data_results = load('results/human6GrossRegions-1-30000.mat');
data_results = load('results/human6AllRegions-1-30000.mat');
human6Results = double(data_results.results);
human6RandomResults = double(data_results.randomResults);
human_gene_info = data_results.gene_info;

data_results = load('results/zapalaMouse-1-30000.mat');
second_dataset_results = data_results.results;
second_dataset_random_results = data_results.randomResults;
second_dataset_genes_info = data_results.gene_info;

addpath('/cortex/code/matlab/homologous_gene_mapping/');
% fprintf('==== MOUSE ====\n');
% gene_to_homolog_group('mouse_laboratory','mouse_laboratory', second_dataset_genes_info.gene_symbols, 'symbol',second_dataset_genes_info.gene_symbols,'symbol',true);
% fprintf('==== HUMAN ====\n' );
% gene_to_homolog_group('human','human', human_gene_info.gene_symbols, 'symbol',human_gene_info.gene_symbols,'symbol',true);
[gene_to_group_matrix_mouse, gene_to_group_matrix_human, homologous_group_id] = gene_to_homolog_group('mouse_laboratory','human', second_dataset_genes_info.gene_symbols, 'symbol',human_gene_info.gene_symbols,'symbol',false);

second_dataset_results = get_group_expression(gene_to_group_matrix_mouse, second_dataset_results');
second_dataset_random_results = get_group_expression(gene_to_group_matrix_mouse, second_dataset_random_results');
human6Results = get_group_expression(gene_to_group_matrix_human, human6Results');
human6RandomResults = get_group_expression(gene_to_group_matrix_human, human6RandomResults');

human_gene_info.gene_symbols = homologous_group_id;
human_gene_info.entrez_ids = homologous_group_id;
second_dataset_genes_info.gene_symbols = homologous_group_id;
second_dataset_genes_info.entrez_ids = homologous_group_id;


data_results = load('results/zapalaMouse_homologs-1-40000.mat');
second_dataset_results = data_results.results;
second_dataset_random_results = data_results.randomResults;
second_dataset_genes_info = data_results.gene_info;



figure('name','Scatter with random');
hold on;
scatterDots = scatterScores(human6Results, human_gene_info, second_dataset_results , second_dataset_genes_info);
randomScatterDots = scatterScores(human6RandomResults, human_gene_info, second_dataset_random_results , second_dataset_genes_info);
set(randomScatterDots, 'CData',random_samples_color);
xlabel('BRO-agreement (ABA6-2013)', 'fontsize',20); ylabel('BRO-agreement Mouse)', 'fontsize',20);
hleg = legend('All','Random');  set(hleg,'Location','Northwest'); legend('boxoff');    set(hleg,'FontSize',20); 
set(gca,'box','off');  
set(gca,'ytick',ytick); set(gca,'xtick',xtick);
title('');
redoScatterImage(xlimit);
saveFigure(gcf, 'figures/scatter_human_mouse_random.png', 'png');
saveFigure(gcf, 'figures/scatter_human_mouse_random.tiff', 'tiff');
saveFigure(gcf, 'figures/scatter_human_mouse_random', 'eps');


% compareDistributions(human6Results, human6RandomResults(:), 'all human6', 'randomHuman6');
% compareDistributions( second_dataset_results,second_dataset_random_results(:), 'all kang','randomKang');
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
ploth = displayMultiDist({second_dataset_results,second_dataset_random_results}, 'BRO-agreement (Kang-2011)',{'All','Random'} ,xlimit);
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

function group_expression = get_group_expression(gene_to_group_matrix, gene_expression)

    num_of_genes_in_group = sum(gene_to_group_matrix,1);
    assert( all(num_of_genes_in_group > 0), 'every group should have atleast one gene');
    
    group_expression = gene_expression * gene_to_group_matrix;
    group_expression = group_expression * diag( 1./ num_of_genes_in_group );
    group_expression = group_expression';
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