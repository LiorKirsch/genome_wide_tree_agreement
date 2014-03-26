function venn_main()

    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

    
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',kangAges);
    geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);       
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',[]);
    load('simpleMeasurementAllRegionsPar.mat');

%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',kangAges);
%     geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);
%     load('simpleMeasurementCortexGrossPar.mat');

    anovaScoresForHuman6 = computeAnova(human_expression, human_gross_region_vec);
    correctedAnova_human6 = mafdr(anovaScoresForHuman6, 'BHFDR', true);
    geneHasFoldChangeLargerThanX_human6 = findMoreThan2Folds(human_expression, human_gross_region_vec);    

    
    anovaScoresForGenesDeveloping = computeAnova(developing_expression, developing_gross_region_vec);
    correctedAnova_developing = mafdr(anovaScoresForGenesDeveloping, 'BHFDR', true);
  


    figure('name','kang venn');
    empirical_pvalue_kang = getEmpiricalPvalues(brainspanResults, brainspanRandomResults);
    empirical_pvalue_kang = mafdr(empirical_pvalue_kang, 'BHFDR', true);
    error = drawVenn(empirical_pvalue_kang < 0.01 , geneHasFoldChangeLargerThanX_dev_all_ages', correctedAnova_developing  < 0.01, 0.5);
    
    figure('name','human6 venn');
    empirical_pvalue_huamn6 = getEmpiricalPvalues(human6Results, human6RandomResults);
    empirical_pvalue_huamn6 = mafdr(empirical_pvalue_huamn6, 'BHFDR', true);
    error = drawVenn(empirical_pvalue_huamn6 < 0.01 , geneHasFoldChangeLargerThanX_human6', correctedAnova_human6  < 0.01, 0.5);
    
    
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