function venn_main()
    init;
    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

%% best resolution cortex regions
    
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',kangAges);
%     geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);       
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',[]);
%     load('data_matfile/simpleMeasurementAllRegionsCortexPar.mat');

%% best resolution regions

    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',kangAges);
    geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);       
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6AllRegions', []);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
%     load('data_matfile/simpleMeasurementAllRegionsPar.mat');
    
    load('simpleMeasurementAllRegions-subjects');
    human6ResultsCell = human6Results;
    human6RandomResultsCell = human6RandomResults;

    num_humans = length(human6ResultsCell);
    for i = 1:num_humans
       
    end



%% 16 gross regions

%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',kangAges);
%     geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);       
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6', []);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
%     load('data_matfile/simpleMeasurement16RegionsPar.mat');
    
%% cortex gross regions

%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',kangAges);
%     geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);
%     load('data_matfile/simpleMeasurementCortexGrossPar.mat');


    num_subjects  = size(human_samples2subjects,2);
    num_genes = size(human_expression,2);
    anovaScoresForHuman6 = nan(num_genes, num_subjects);
    correctedAnova_human6 = nan(num_genes, num_subjects);
    meanHuman6Results = zeros(size(human6ResultsCell{1}) );
    meanHuman6RandomResults = [];
  
    for i = 1:num_subjects
        current_subject_samples = human_samples2subjects(:,i);
        anovaScoresForHuman6(:,i) = computeAnova(human_expression(current_subject_samples,:) , human_gross_region_vec(current_subject_samples) );
        correctedAnova_human6(:,i) = mafdr( anovaScoresForHuman6(:,i), 'BHFDR', true);
        meanHuman6Results = meanHuman6Results + human6ResultsCell{i} / num_humans;
        meanHuman6RandomResults = [meanHuman6RandomResults , human6RandomResultsCell{i}] ;
     end
    
%     anovaScoresForHuman6 = computeAnova(human_expression, human_gross_region_vec);
%     correctedAnova_human6 = mafdr(anovaScoresForHuman6, 'BHFDR', true);
%     geneHasFoldChangeLargerThanX_human6 = findMoreThan2Folds(human_expression, human_gross_region_vec);    

    
%     anovaScoresForGenesDeveloping = computeAnova(developing_expression, developing_gross_region_vec);
%     correctedAnova_developing = mafdr(anovaScoresForGenesDeveloping, 'BHFDR', true);
  


%     figure('name','human6 venn');
    empirical_pvalue_huamn6 = getEmpiricalPvalues(meanHuman6Results, meanHuman6RandomResults);
    empirical_pvalue_huamn6 = mafdr(empirical_pvalue_huamn6, 'BHFDR', true);
    BRO_significant = empirical_pvalue_huamn6 < 0.01;
    ANOVA_significant = mean(correctedAnova_human6,2)  < 0.01;
    fprintf('BRO but not ANOVA %d, (%g)\n',  sum(BRO_significant) - sum(BRO_significant& ANOVA_significant), (sum(BRO_significant) - sum(BRO_significant& ANOVA_significant))/ num_genes);
    printGeneNames( (BRO_significant & ~ANOVA_significant) , human_gene_info.gene_symbols ,'BRO_not_ANOVA.txt');
    fprintf('ANOVA but not BRO %d, (%g)\n',  sum(ANOVA_significant) - sum(BRO_significant& ANOVA_significant), (sum(ANOVA_significant) - sum(BRO_significant& ANOVA_significant))/ num_genes);
    fprintf('ANOVA and BRO %d, (%g)\n',  sum(BRO_significant& ANOVA_significant), sum(BRO_significant& ANOVA_significant)/ num_genes  );
%     error = drawVenn(empirical_pvalue_huamn6 < 0.01 , correctedAnova_human6  < 0.01, 0.5);
    error = vennX( [ sum(BRO_significant), sum(BRO_significant& ANOVA_significant), sum(ANOVA_significant) ], 0.5 );
    
    
%     figure('name','kang venn');
%     empirical_pvalue_kang = getEmpiricalPvalues(brainspanResults, brainspanRandomResults);
%     empirical_pvalue_kang = mafdr(empirical_pvalue_kang, 'BHFDR', true);
%     error = drawVenn(empirical_pvalue_kang < 0.01 , geneHasFoldChangeLargerThanX_dev_all_ages', correctedAnova_developing  < 0.01, 0.5);
    
   
    
end

function printGeneNames(booleanVector, geneSym,fileName)
    fid = fopen(fileName, 'w');
    geneSymSubset = geneSym(booleanVector);
    for i =1:length(geneSymSubset)
        fprintf(fid,'%s\n', geneSymSubset{i} );
    end
    fclose(fid);
end
function printGeneNamesOrdered(geneScore, geneSym,fileName)
% direction should be 'descend' or 'ascend'

    [sortedScore, sortInd] = sort(geneScore, direction);
    fid = fopen(fileName, 'w');
    sortedGeneSym = geneSym(sortInd);
    for i =1:length(sortedGeneSym)
        fprintf(fid,'%g,%s\n', sortedScore(i), sortedGeneSym{i} );
    end
    fclose(fid);
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