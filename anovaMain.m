function anovaMain(startIndex, finishIndex)

    addpath('/home/lab/gal/develop/matlab');
    addpath('~/Projects/general use functions/');

    % Load the tree
    data_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability');
    filename = 'humanOntologyObject.mat';
    load(fullfile(data_dirname, filename), 'humanOntology'); %#ok

    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',kangAges);    
    geneHasFoldChangeLargerThanX_dev_all_ages = findMoreThan2Folds(developing_expression, developing_gross_region_vec);    
    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6AllRegions',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);
        [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
%     grossStructures = {'Frontal Lobe','hippocampal formation', 'Occipital Lobe', 'Parietal Lobe','Temporal Lobe','Amygdala','Striatum','Thalamus','Cerebellum'};

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',[]);

    %normalize the data to have zero mean and unit variance
%     human_expression = human_expression - repmat(mean(human_expression,2),1,size(human_expression,2) );
%     human_expression = human_expression ./ repmat( std(human_expression,0,2),1,size(human_expression,2) );
%     developing_expression = developing_expression - repmat(mean(developing_expression,2),1,size(developing_expression,2) );
%     developing_expression = developing_expression ./ repmat( std(developing_expression,0,2),1,size(developing_expression,2) );
    %-----------------------------------------------------%
    
    num_subjects  = size(human_samples2subjects,2);
    num_genes = size(human_expression,2);
    anovaScoresForHuman6 = nan(num_genes, num_subjects);
    correctedAnova_human6 = nan(num_genes, num_subjects);
    for i = 1:num_subjects
        current_subject_samples = human_samples2subjects(:,i);
        anovaScoresForHuman6(:,i) = computeAnova(human_expression(current_subject_samples,:) , human_gross_region_vec(current_subject_samples) );
        correctedAnova_human6(:,i) = mafdr( anovaScoresForHuman6(:,i), 'BHFDR', true);
    end
    geneHasFoldChangeLargerThanX_human6 = findMoreThan2Folds(human_expression, human_gross_region_vec);    
    
    anovaScoresForGenesDeveloping = computeAnova(developing_expression, developing_gross_region_vec);
    correctedAnova_developing = mafdr(anovaScoresForGenesDeveloping, 'BHFDR', true);
    
%     load('kang_genes_2_fold.mat','geneHasFoldChangeLargerThanX');
%     save('kang_genes_2_fold.mat','geneHasFoldChangeLargerThanX_dev_all_ages','correctedAnova_developing','geneHasFoldChangeLargerThanX_human6','correctedAnova_human6');
    anovaScoresForGenesDeveloping_withFold = computeAnova(developing_expression(:,geneHasFoldChangeLargerThanX_dev_all_ages), developing_gross_region_vec);
    correctedAnova_developing_withFold = mafdr(anovaScoresForGenesDeveloping_withFold, 'BHFDR', true);    
    

%     anovaScoresForGenesHuman6 = computeAnova(human_expression, human_gross_region_vec);
%     correctedAnova_human6 = mafdr(anovaScoresForGenesHuman6, 'BHFDR', true);

    save('anova_results_cortex.mat', 'correctedAnova_developing_withFold', 'correctedAnova_developing','geneDexLarger2','developing_genes_info');
end


function anovaScoresForGenes = computeAnova(expressionMatrix, groupData)
    numOfGenes = size(expressionMatrix,2);
    num_of_samples = size(expressionMatrix,1);
    anovaScoresForGenes = nan(numOfGenes,1);
    parfor i =1:numOfGenes %number of genes
        anovaScoresForGenes(i) = anova1(expressionMatrix(:,i), groupData,'off') ;
    end

end

function treeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology)
%     treeDistances = humanOntology.meanDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
    treeDistances = humanOntology.unDirectedDistanceMatrix(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.longDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.shortDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
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

function samplesDistance = getSamplesDistance(sample_gross_ind, distance_matrix)

    numberOfSamples = length(sample_gross_ind);
    samplesDistance = nan(numberOfSamples,numberOfSamples);
    
    for  i =1:numberOfSamples
        for  j =1:numberOfSamples
           area_sample_i = sample_gross_ind(i);
           area_sample_j = sample_gross_ind(j);
           samplesDistance(i,j) =  distance_matrix(area_sample_i, area_sample_j);
        end
    end
           
end


