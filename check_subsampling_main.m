function check_subsampling_main()

    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

    
%     load('data_matfile/simpleMeasurement16RegionsPar.mat');
    load('data_matfile/simpleMeasurementAllRegionsPar.mat');

     

    figure('name','kang venn');
    empirical_pvalue_kang = getEmpiricalPvalues(brainspanResults, brainspanRandomResults);
    empirical_pvalue_kang = mafdr(empirical_pvalue_kang, 'BHFDR', true);
    
    empirical_pvalue_huamn6 = getEmpiricalPvalues(human6Results, human6RandomResults);
    empirical_pvalue_huamn6 = mafdr(empirical_pvalue_huamn6, 'BHFDR', true);

    
    num_repeat = 10;
    empirical_pvalue_huamn6_subset = nan(length(empirical_pvalue_huamn6), num_repeat);
    for i =1:num_repeat
        BRO_data = load(['results/simpleMeasurementAllRegions-1-3000',num2str(i-1),'.mat']');
        empirical_pvalue_huamn6_subset(:,i) = getEmpiricalPvalues(BRO_data.human6Results, BRO_data.human6RandomResults);
        empirical_pvalue_huamn6_subset(:,i) = mafdr(empirical_pvalue_huamn6_subset(:,i), 'BHFDR', true);
    end
    
    figure('name','BRO scores');
    plot_bars(empirical_pvalue_huamn6, empirical_pvalue_huamn6_subset, empirical_pvalue_kang);
end

function plot_bars(empirical_pvalue_huamn6, empirical_pvalue_huamn6_subset, empirical_pvalue_kang)

    
    y = [sum(empirical_pvalue_huamn6 < 0.01),  mean(sum(empirical_pvalue_huamn6_subset < 0.01,1)),    sum(empirical_pvalue_kang < 0.01)];
    
    e = [0,  std(sum(empirical_pvalue_huamn6_subset < 0.01,1)),    0];
    hold on;
    bar(y);
    errorbar(y,e,'xr')
    hold off;
    xticklabel_rotate(1:3,45,{'human6','human6 equal size','kang'});
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

function geneHasFoldChangeLargerThanX = findMoreThan2Folds(geneExpression, region_vector)
    fold_change_threshold = 2;
    numberOfGenes = size(geneExpression,2);
    numberOfSamples = size(geneExpression,1);
    
    assert(length(region_vector) == numberOfSamples);
    
    numberOfRegions = length( unique(region_vector) );
    
    meanExpressionInRegion = nan(numberOfRegions, numberOfGenes);
    for i = 1:numberOfRegions
        meanExpressionInRegion(i,:) = mean(geneExpression(region_vector == i,:),1);
    end
    
    
    geneHasFoldChangeLargerThanX = false(1,numberOfGenes);
    for i = 1:numberOfRegions
        for j = i+1:numberOfRegions
            fold_change = abs( meanExpressionInRegion(i,:) - meanExpressionInRegion(j,:) );
            geneHasFoldChangeLargerThanX = geneHasFoldChangeLargerThanX | (fold_change > log(fold_change_threshold));
        end
    end

end


function anovaScoresForGenes = computeAnova(expressionMatrix, groupData)
    numOfGenes = size(expressionMatrix,2);
    num_of_samples = size(expressionMatrix,1);
    anovaScoresForGenes = nan(numOfGenes,1);
    parfor i =1:numOfGenes %number of genes
        anovaScoresForGenes(i) = anova1(expressionMatrix(:,i), groupData,'off') ;
    end

end