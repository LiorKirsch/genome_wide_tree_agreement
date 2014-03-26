function treeOrderingMainSimpleKangAge(startIndex, finishIndex)

    addpath('/home/lab/gal/develop/matlab');

    % Load the tree
    addpath('/home/lab/lior/Projects/buildStructureOntology/');
    filename = '/home/lab/lior/Projects/buildStructureOntology/humanOntologyObjectV2.mat';
    load(filename, 'humanOntology'); %#ok

    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

    for i =4:length(kangAges)
        [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',kangAges{i});

% %         normalize the data to have zero mean and unit variance
%         developing_expression = developing_expression - repmat(mean(developing_expression,2),1,size(developing_expression,2) );
%         developing_expression = developing_expression ./ repmat( std(developing_expression,0,2),1,size(developing_expression,2) );
% %         -----------------------------------------------------%
        if any(strcmp(developing_gross_region_names, 'Hippocampus'))
            developing_gross_region_names{ strcmp(developing_gross_region_names, 'Hippocampus') } = 'hippocampal formation';
        end

        [~,gross_developing_region_indices_in_ontology] = ismember(developing_gross_region_names, humanOntology.structureLabels(:,4) );
        developingTreeDistances = getTreeDistance(humanOntology, gross_developing_region_indices_in_ontology);
        samplesDevelopingTreeDistance = getSamplesDistance(developing_gross_region_vec, developingTreeDistances);

        [brainspanResults{i},brainspanRandomResults{i}] = compareExpressionToDistanceMetric(developing_expression, samplesDevelopingTreeDistance ,1, 30000);
        fprintf('\n-%d-\n',i);

    end
    fileName = 'results/simpleMeasurementKangAges.mat';
    save(fileName, 'brainspanResults','brainspanRandomResults','developing_genes_info');
end



function drawDistribution(brainspanResults, kangAges)
    dataMerged = nan(length(brainspanResults{4}), length(brainspanResults) - 3);
    for i =4:length(brainspanResults)
        dataMerged(:,i) = brainspanResults{i};
    end
    [dataInBins ,xout] = hist(dataMerged,100);
    plot(xout,dataInBins);
    legend(kangAges{4:end});
end

function treeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology)
%     treeDistances = humanOntology.meanDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
    treeDistances = humanOntology.unDirectedDistanceMatrix(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.longDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.shortDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
end

function [result, random_result] = compareExpressionToDistanceMetric( human_expression, distancesMatrix ,startIndex, finishIndex)
    

    
    [numberOfSamples, numberOfGenes] = size(human_expression);
    
    rand('state', 1221);
    
    repeat_random = 20;
    randomPerm = nan(numberOfSamples , repeat_random);
    for i = 1:repeat_random
        randomPerm(:,i) = randperm(numberOfSamples);
    end
    random_human_expression = human_expression(randperm(numberOfSamples),:);
    
    result = nan(numberOfGenes,1);
    random_result = nan(numberOfGenes,repeat_random);
    
    fprintf('progress:    ');

    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distancesMatrix(onlyUpperTri);
    
    onlyUpperRandomDistanceMatrix =  nan(length(onlyUpperDistanceMatrix) ,repeat_random);
    for i = 1:repeat_random
        randomDistanceMatrix = distancesMatrix( randomPerm(:,i), randomPerm(:,i) );
        onlyUpperRandomDistanceMatrix(:,i) = randomDistanceMatrix(onlyUpperTri);
    end
    
    lastIndex = min(finishIndex, numberOfGenes);
    for i = startIndex: lastIndex
        current_gene_expression = human_expression(:,i);
        expression_distance_matrix = squareform( pdist(current_gene_expression,'euclidean') );
        onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
        result(i) = corr(onlyUpperExpressionMatrix, onlyUpperDistanceMatrix , 'type','Spearman');
%         random_result(i,:) = corr(onlyUpperExpressionMatrix, onlyUpperRandomDistanceMatrix , 'type','Spearman');
        printPercentCounter(i,numberOfGenes);
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

