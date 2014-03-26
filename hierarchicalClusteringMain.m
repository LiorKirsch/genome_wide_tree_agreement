function hierarchicalClusteringMain()

%     conf = configureation();
    
    %============== Load the tree ==================
    addpath('~/Projects/individual variability');
    addpath('~/Projects/general use functions/');
    addpath('~/Projects/genome_wide_agreement/');
    addpath('/home/lab/gal/develop/matlabxt');

    load('~/Projects/individual variability/humanOntologyObject.mat');
    treeMatrix = humanOntology.dependencyMatrix;
    
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6', []);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);

    region_reorder = [1     7     2     5     6     4     3     8    10     9    11    13    15    16    12    14];
    [correlationMatrixReorder,samples_region_reorder, samples2subjects] = createMeanPairCorrelationMatrixWithOrder(human_expression, human_gross_region_vec, region_reorder, human_samples2subjects);
    
    drawCorrelationMatrixWithBorders(correlationMatrixReorder, samples_region_reorder, human_gross_region_names(region_reorder));
    colorbar;
    axis ij;
    
%     correlationMatrix = createMeanPaircorrelationMatrix(human_expression, human_gross_region_vec, length( human_gross_region_names));
%     
%     distanceMatrix = 1 - correlationMatrix ;
%     upperHalf = isnan(distanceMatrix);
%     distanceMatrix(upperHalf) = 0;
%     distanceMatrix = distanceMatrix + distanceMatrix';
%     y = squareform(distanceMatrix);
%     links  = linkage(y);
%     figure;
%     subplot(2,1,2);
%     [~,~,outperm] = dendrogram(links );
%     xticklabel_rotate(1:length(human_gross_region_names), 40, human_gross_region_names(outperm));
%     
%     subplot(2,1,1);
%     permutedDistances = distanceMatrix(outperm,outperm);
%     permutedDistances( upperHalf) = nan;
%     imagesc(1 - permutedDistances);
%     xticklabel_rotate(1:length(human_gross_region_names), 40, human_gross_region_names(outperm));
%     colorbar()
    
    mean_expression_in_region = nan(length( human_gross_region_names), size(human_expression,2));
    median_expression_in_region = nan(length( human_gross_region_names), size(human_expression,2));
    for i = 1:length(human_gross_region_names)
        mean_expression_in_region(i,:) = mean( human_expression( human_gross_region_vec == i, :),1);
        median_expression_in_region(i,:) = median( human_expression( human_gross_region_vec == i, :),1);
    end
    
    createFigure();
    a = linkage(mean_expression_in_region, 'average');
    [~,~,outperm] = dendrogram(a);
    xticklabel_rotate(1:length(human_gross_region_names), 40, human_gross_region_names(outperm));
    saveFigure(gcf, 'hierarachical-mean', 'png');
    
    createFigure();
    a = linkage(median_expression_in_region);
    [~,~,outperm] = dendrogram(a);
    xticklabel_rotate(1:length(human_gross_region_names), 40, human_gross_region_names(outperm));
    saveFigure(gcf, 'hierarachical-media', 'png');
end

function [correlationMatrix,samples_region, samples2subjects] = createMeanPairCorrelationMatrixWithOrder(human_expression, human_region_vec, region_order, samples2subjects)

    
    all_pairs_corr = corr(human_expression');
    [~, region_order_ind] = sort(region_order);
    new_sample_order = region_order_ind(human_region_vec);
    [~, sort_ind] = sort(new_sample_order);
    
    correlationMatrix = all_pairs_corr(sort_ind,sort_ind);
    samples_region = human_region_vec(sort_ind);
    samples2subjects = samples2subjects(sort_ind,:);
end

function drawCorrelationMatrixWithBorders(correlationMatrix, sample_group,labels)

    number_of_samples = length(sample_group);
    borders = ~(diff(sample_group) == 0) ;
    borders = find ([false; borders]);
    hold on;
    imagesc(correlationMatrix);
%     for i = 1: length(borders)
%         borderIndex = borders(i);
%         plot(1:number_of_samples, ones(1,number_of_samples) * borderIndex,'w');
%         plot(ones(1,number_of_samples) * borderIndex, 1:number_of_samples,'w');
%     end
    axis equal;
    xticklabel_rotate([borders; number_of_samples], 40, labels);

    hold off;
    hold off;
    
end

