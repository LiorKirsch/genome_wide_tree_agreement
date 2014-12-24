function viewScatterPlotsSubjects()


%     load('simpleMeasurement16Regions-subjects');

load('/home/lab/lior/Projects/genome_wide_agreement/results/simpleMeasurementAllRegions-1-30000-subjects.mat')

%     load('simpleMeasurementAllRegions-subjects');
%     load('simpleMeasurementAllCortexRegions-subjects.mat');
    human6ResultsCell = human6Results;
    human6RandomResultsCell = human6RandomResults;
    
    
    drawAllSubjectsInOnePlot(human6ResultsCell,  cell2mat(human6RandomResultsCell' ));
    
    
    countGenesAboveThreshAllSubjects(human6ResultsCell, human6RandomResultsCell);
    countGenesAboveThreshAllSubjects( brainspanResults , brainspanRandomResults);

    num_humans = length(human6ResultsCell);
     
    Human6ResultsTogether = cell2mat(human6ResultsCell');
    meanHuman6RandomResults = cell2mat(human6RandomResultsCell');
    meanHuman6Results = mean(Human6ResultsTogether,2);
    stdHuman6Results = std(Human6ResultsTogether,1,2);
    
   
    createAllFigures(meanHuman6Results, meanHuman6RandomResults, human_gene_info, 0 );
    
   

    for i = 1:num_humans
        createAllFigures(human6ResultsCell{i}, human6RandomResultsCell{i}, human_gene_info,i);
    end


end

function drawAllSubjectsInOnePlot(human6ResultsCell, human6RandomResults)
    % xlimit = [0.4,0.9];
    xlimit = [-0.3,0.8]; % full tree
%     xlimit = [-0.1,0.25]; % for cortex
    ylimit = [0, 0.25] ; % full tree
%     ylimit = [0, 0.4] ; % for cortex

color = nan(2,3);
random_samples_color = [0.4, 0.4,0.4];
color(1,:) = random_samples_color;
color(2,:) =  [0.8, 0.8,0.8];
color(3,:) =  [0.4,0.7,0] ;
color(4,:) =  [ 0.7, 0.4, 0.9];
color(5,:) =   [ 0.9, 0.2, 0.5];

% colors = distinguishable_colors(20);
colors = [color; [ 0 0 1;  1 0 0 ; 0.5,0.8,0.3;0 0 0.3;] ];

    folder = 'seperateSubjects';

    figure('Name' ,'samples vs random');
    desc = {'Random','Human1', 'Human2', 'Human3', 'Human4', 'Human5', 'Human6'};
    ploth = displayMultiDist( [human6RandomResults(:); human6ResultsCell ],'BRO-agreement (ABA6-2013)', desc ,xlimit);
    for i = 1:length(desc)
        set(ploth(i), 'Color', colors(i,:) );
        set(ploth,'LineWidth',  3 );
        set(ploth,'LineStyle',  '-' );
    end
    set(ploth(1),'LineStyle',  '-.' );

    
    redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
    set(gca,'ytick',[0:0.2:1]);
    ylim( [0, 0.6]);
    
    fileName = fullfile(folder,'joint_figure_dist');
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');
    
    
end

function createAllFigures(human6Results, human6RandomResults, human_gene_info, index)

folder = 'seperateSubjects';
random_samples_color = [0.4, 0.4,0.4];
all_samples_color = [0.8, 0.8,0.8];
house_keeping_color = [0.4,0.7,0] ;
hox_color = [ 0.7, 0.4, 0.9];
axon_guidance_color =  [ 0.9, 0.2, 0.5];

% colors = distinguishable_colors(20);
colors = [all_samples_color; [ 0 0 1;  0 0 0.3; 1 0 0 ; 0.5,0.8,0.3] ];

set(0,'DefaultAxesColorOrder',colors);
set(0,'DefaultAxesLineStyleOrder',{'-',':','--','.-'});

xtick = -0.1:0.1:1;
ytick = -0.1:0.1:1;


    figure('Name','qvalue vs gene count');
    ploth = drawNumGenesVsRightPvalue(human6Results,human6RandomResults);
%     set(gca,'ytick',[0:5000:20000]);
    fileName = fullfile(folder,sprintf('pvalue_cumsumm_%d',index ) );
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');
    

    %=============== counts how many genes are below 0.01 q-value============
    % empirical_pvalue_human6 = getEmpiricalPvalues(human6Results, human6RandomResults);
    % empirical_pvalue_human6 = mafdr(empirical_pvalue_human6, 'BHFDR', true);
    % fprintf('all genes: %d / %d  (%d%%)\n', sum(empirical_pvalue_human6 < 0.01) , length(empirical_pvalue_human6) ,floor(sum(empirical_pvalue_human6< 0.01) / length(empirical_pvalue_human6)*100) );
    % empirical_pvalue_developing = getEmpiricalPvalues(brainspanResults, brainspanRandomResults);
    % empirical_pvalue_developing = mafdr(empirical_pvalue_developing, 'BHFDR', true);
    % fprintf('developing: %d / %d  (%d%%)\n', sum(empirical_pvalue_developing < 0.01) , length(empirical_pvalue_developing) ,floor(sum(empirical_pvalue_developing< 0.01) / length(empirical_pvalue_developing)*100) );
    %========================================================================

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

    moreThenXPercent(human6Results, human6RandomResults, 0.99);

    paxGeneSymbol = {'Pax1','Pax2','Pax3','Pax4','Pax5','Pax6','Pax7','Pax8';};
    [~, indInList] = ismember(upper(paxGeneSymbol),human_gene_info.gene_symbols);
    paxGene.symbol = human_gene_info.gene_symbols(indInList);
    paxGene.entrez = human_gene_info.entrez_ids(indInList);
    [paxGene.scores, paxGene.symbols, paxGene.entrez] = addScoresToSubset(human6Results, human_gene_info, paxGene.entrez, paxGene.symbol);
  

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
    compareDistributions(housekeeping_subset.scores, human6Results,  'house', 'all human6');

    % compareDistributions( hoxGenes_subset.scores, human6RandomResults(:), 'hox','randomHuman6');
    % compareDistributions( housekeeping_subset.scores, human6RandomResults(:), 'house','randomHuman6');

    % xlimit = [0.4,0.9];
    xlimit = [-0.3,0.8]; % full tree
%     xlimit = [-0.1,0.25]; % for cortex
    ylimit = [0, 0.25] ; % full tree
%     ylimit = [0, 0.4] ; % for cortex

    f = figure('name','Human6 distribution');  
    ploth = displayMultiDist({human6Results,human6RandomResults}, 'BRO-agreement (ABA6-2013)',{'All','Random'} ,xlimit);
    set(ploth(2), 'Color', random_samples_color );  set(ploth,'LineWidth',  3 );
    redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
    set(gca,'ytick',[0:0.2:1]);
    ylim( [0, 0.6]);
    fileName = fullfile(folder,sprintf('dist_random_%d',index ) );
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');

   
    figure('name','Human6 cahoy distribution'); 
    ploth = displayMultiDist({human6Results,cahoy_human_subset.scores_neurons,  cahoy_human_subset.scores_oligo, cahoy_human_subset.scores_astro}, 'BRO-agreement (ABA6-2013)',{'All','Neurons','Oligodendrocytes','Astrocytes'} ,xlimit);
    redoDistImage(xlimit); set(ploth,'LineWidth',  3 );
    ylim(ylimit);
    fileName = fullfile(folder,sprintf('dist_cellType_%d',index) );
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');

    figure('name','Human6 housekeeping distribution'); 
    ploth = displayMultiDist({human6Results,housekeeping_subset.scores}, 'BRO-agreement (ABA6-2013)',{'All','Housekeeping'} ,xlimit);
    set(ploth(2),'Color', house_keeping_color );  set(ploth,'LineWidth',  3 );
    redoDistImage(xlimit);
    ylim(ylimit);
    fileName = fullfile(folder,sprintf('dist_house_%d',index) );
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');
    
    figure('name','Human6 axon guidance distribution'); 
    ploth = displayMultiDist({human6Results,axon_guidance_subset.scores,hoxGenes_subset.scores,paxGene.scores}, 'BRO-agreement (ABA6-2013)',{'All','Axon guidance','Hox','PAX6'} ,xlimit);
    set(ploth(2),'Color', axon_guidance_color ); 
    set(ploth(3),'Color', hox_color ); 
    set(ploth,'LineWidth',  3 );
    redoDistImage(xlimit);
    ylim(ylimit);
    fileName = fullfile(folder,sprintf('dist_axon_%d',index) );
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'tiff');
    saveFigure(gcf, fileName, 'eps');
end

function countGenesAboveThreshAllSubjects(human6ResultsCell, human6RandomResultsCell)
    num_genes = size(human6ResultsCell{1},1);
    num_subjects = length(human6ResultsCell);

    geneSig = nan(num_genes, num_subjects);

    allRandomResults = [];
    allResults = [];
    for i = 1:num_subjects
        human6Results = human6ResultsCell{i};
        allResults = [allResults, human6Results];
        human6RandomResults = human6RandomResultsCell{i};
        allRandomResults = [allRandomResults, human6RandomResults];
        fprintf('subject %d: ',i);
        geneLargerThanX = moreThenXPercent(human6Results, human6RandomResults, 0.99);
        geneSig(:,i) = geneLargerThanX;
    end

    pairs = nchoosek(1:num_subjects,2); pairs = [pairs;  [1:num_subjects; 1:num_subjects]'];
    correMatrix = nan(num_subjects);
    jointMatrix = nan(num_subjects);
    for i=1:size(pairs,1)
       correlationPair = corr( allResults(:, pairs(i,1)), allResults(:, pairs(i,2) ),'type','Spearman');
       correMatrix( pairs(i,1) , pairs(i,2) ) = correlationPair;
       correMatrix( pairs(i,2) , pairs(i,1) ) = correlationPair;
       
       jointMatrix( pairs(i,1) , pairs(i,2) ) = sum( geneSig(:,pairs(i,1)) & geneSig(:,pairs(i,2)) ) / num_genes;
       jointMatrix( pairs(i,2) , pairs(i,1) ) = sum( geneSig(:,pairs(i,1)) & geneSig(:,pairs(i,2)) ) / num_genes;
    end
%     correMatrix( logical(eye(size(correMatrix))) ) = nan;
    figure('name','Correlations between subject'); 
    set(gca,'Fontsize',20);
    h = imagesc(correMatrix); colorbar;
%     set(h,'alphadata',~isnan(correMatrix))

    xlabel('subject index','fontsize',20); ylabel('subject index','fontsize',20);
    saveFigure(gcf, sprintf('subject_score_correlation.png'), 'png');
    saveFigure(gcf, sprintf('subject_score_correlation.tiff'), 'tiff');
    correlationPairs = correMatrix( triu(true(size(correMatrix) ),1 ));
    fprintf('correlation across pairs of subjects %g, (+- %g)\n',mean( correlationPairs(:) ) , std( correlationPairs(:) ) );
    figure('name','Number of joint genes above threshold'); 
    imagesc(jointMatrix); colorbar;
    set(gca,'Fontsize',20);
    xlabel('subject index','fontsize',20); ylabel('subject index','fontsize',20);
    how_many_larger = sum(  all(geneSig,2) );
    fprintf( 'All subjects - %d%% (%d) of the genes are larger \n', round( how_many_larger / num_genes*100), how_many_larger);

    fprintf('using mean - ');
    geneLargerThanX = moreThenXPercent(mean(allResults,2), allRandomResults, 0.99);
    fprintf('using median - ');
    geneLargerThanX = moreThenXPercent(median(allResults,2), allRandomResults, 0.99);
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

function scatterDots = scatterScores(scoresA, geneInfoA, scoresB, geneInfoB)
    [ind_for_A, ind_for_B] = mapGenes(geneInfoA.gene_symbols, geneInfoB.gene_symbols, geneInfoA.entrez_ids, geneInfoB.entrez_ids);
    intersect_scoresA = scoresA(ind_for_A);
    intersect_scoresB = scoresB(ind_for_B);
    scatterDots = scatter(intersect_scoresA, intersect_scoresB,1,'filled');
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
    xlabel('q-value (FDR)','Fontsize',20);
    ylabel('number of genes','Fontsize',20);
    set(gca,'Fontsize',20); 
    xlim([-0.05 1.0]);
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

function geneLargerThanX = moreThenXPercent(scores, randomScores, x_percent)

% randomScoreSorted = sort(randomScores(:) );
% 
% index = round( x_percent * length(randomScoreSorted) );
% geneLargerThanX = scores > randomScoreSorted(index);
% how_many_larger = sum(  geneLargerThanX );

x_percent = 1  - x_percent;
empirical_pvalue = getEmpiricalPvalues(scores, randomScores(:) );
geneLargerThanX = empirical_pvalue < x_percent;
how_many_larger = sum(  geneLargerThanX );

fprintf( '%d%% (%d) of the genes are larger \n', round( how_many_larger / length(scores)*100), how_many_larger);
end

function redoDistImage(xlimit)
    xtick = -0.2:0.2:1;
    ytick = -0.2:0.1:1;

    set(gca,'Fontsize',20);         %set(gca,'xscale','log');
    set(gca,'ytick',ytick); set(gca,'xtick',xtick);
    xlim( xlimit)
end

function redoScatterImage()
    xtick = -0.2:0.2:1;
    ytick = -0.2:0.2:1;

    set(gca,'Fontsize',20);         %set(gca,'xscale','log');
    set(gca,'ytick',ytick); set(gca,'xtick',xtick);
    xlim( [-0.13,0.8])
    ylim( [-0.13,0.8])
    axis equal
end

function scatterDots = scatterSubsetScores(scoresA, geneInfoA, scoresB, geneInfoB, subset_entrez, subset_names,colorTriplet)
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
%      set(scatterDots,'FaceColor',colorTriplet);
    
    
end

function [scoresA, gene_symbols, entrez_ids] = addScoresToSubset(scoresA, geneInfo, subset_entrez, subset_names)
    [ind_for_A, ~] = mapGenes(geneInfo.gene_symbols,subset_names, geneInfo.entrez_ids, subset_entrez);
        
    scoresA = scoresA(ind_for_A);
    gene_symbols = geneInfo.gene_symbols(ind_for_A);
    entrez_ids = geneInfo.entrez_ids(ind_for_A);
end

function ploth = displayMultiDist(scores_cell, xlabel_string,legendStrings,range)

    num_bins = 50;

    spacing = linspace(range(1),range(2),num_bins);
    vals_in_bins = nan(length(spacing), length(scores_cell));
    for i =1 :length(scores_cell);
        scores = scores_cell{i};
        vals_in_bins(:,i) = histc(scores(:), spacing);
%         vals_in_bins(:,i) = smooth(vals_in_bins(:,i), 0.1, 'moving');
        vals_in_bins(:,i) = vals_in_bins(:,i) / sum(vals_in_bins(:,i));
    end

    vals_in_bins_smooth = smooth(vals_in_bins,0.1 / num_bins);
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