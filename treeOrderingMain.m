function treeOrderingMain()

    addpath('/home/lab/gal/develop/matlab');

    % Load the tree
    data_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability');
    filename = 'humanOntologyObject.mat';
    load(fullfile(data_dirname, filename), 'humanOntology'); %#ok

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6AllRegions',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);
%     grossStructures = {'Frontal Lobe','hippocampal formation', 'Occipital Lobe', 'Parietal Lobe','Temporal Lobe','Amygdala','Striatum','Thalamus','Cerebellum'};

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);

    if any(strcmp(developing_gross_region_names, 'Hippocampus'))
        developing_gross_region_names{ strcmp(developing_gross_region_names, 'Hippocampus') } = 'hippocampal formation';
    end
    [~,gross_region_indices_in_ontology] = ismember(human_gross_region_names, humanOntology.structureLabels(:,4) );
    [~,gross_developing_region_indices_in_ontology] = ismember(developing_gross_region_names, humanOntology.structureLabels(:,4) );
    
    human6TreeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology);
    brainspanTreeDistances = getTreeDistance(humanOntology, gross_developing_region_indices_in_ontology);
    
    human6FlatRegionDistances = squareform( pdist((1:length(gross_region_indices_in_ontology) )') );
    brainspanFlatRegionDistances = squareform( pdist((1:length(gross_developing_region_indices_in_ontology) )') );
    human6PhysicalDistances = squareform( pdist(physicalLocation,'euclidean') );
    
    
    identity_region_mapping = (1:length(human_gross_region_vec))';
    
    [human6Results,human6RandomResults] = compareExpressionToDistanceMetric(human_expression , human6TreeDistances, human_gross_region_vec);
    [brainspanResults,brainspanRandomResults] = compareExpressionToDistanceMetric(developing_expression, brainspanTreeDistances , developing_gross_region_vec);
    
    [human6FlatResults,human6FlatRandomResults] = compareExpressionToDistanceMetric(human_expression , human6FlatRegionDistances, human_gross_region_vec);
    [brainspanFlatResults,brainspanFlatRandomResults] = compareExpressionToDistanceMetric(developing_expression, brainspanFlatRegionDistances , developing_gross_region_vec);

    [human6PhysicalResults,human6PhysicalRandomResults] = compareExpressionToDistanceMetric(human_expression , human6PhysicalDistances, identity_region_mapping);

    save('cortexAllTreeResultsTriplets.mat', 'human6Results', 'human6RandomResults','human_gene_info','brainspanResults','brainspanRandomResults','developing_genes_info', 'brainspanFlatResults','brainspanFlatRandomResults' ,'human6FlatResults','human6FlatRandomResults' ,'human6PhysicalResults','human6PhysicalRandomResults');
end

function treeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology)

%     treeDistances = humanOntology.meanDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
    treeDistances = humanOntology.unDirectedDistanceMatrix(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.longDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = humanOntology.shortDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);

end

function [result, random_result] = compareExpressionToDistanceMetric( human_expression, distancesMatrix, gross_region_vec )
    

    
    [numberOfSamples, numberOfGenes] = size(human_expression);
    tripletsIndices = createRandomTriplets(numberOfSamples,10^6); 
    
    rand('state', 1221);
    
    repeat_random = 100;
    randomPerm = nan(numberOfSamples , repeat_random);
    for i = 1:repeat_random
        randomPerm(:,i) = randperm(numberOfSamples);
    end
    random_human_expression = human_expression(randperm(numberOfSamples),:);
    
    result = nan(numberOfGenes,1);
    random_result = nan(numberOfGenes,repeat_random);
    
    [treeDistance12, treeDistance13] = calcTripletTreeDistance(distancesMatrix,gross_region_vec, tripletsIndices);
    
    treeDistance12_random = nan(length(treeDistance12),repeat_random);
    treeDistance13_random = nan(length(treeDistance13),repeat_random);
    fprintf('permuted triplet distance generation:    ');
    for i = 1:repeat_random
        [treeDistance12_random(:,i), treeDistance13_random(:,i)] = calcTripletTreeDistance(distancesMatrix,gross_region_vec(randomPerm(:,i)), tripletsIndices);
        printPercentCounter(i,repeat_random);
    end
    fprintf('\n');
    
    distance12BiggerInTree = treeDistance12 > treeDistance13;
    distance12BiggerInTree_random = treeDistance12_random > treeDistance13_random;
    fprintf('progress:    ');

    
    parfor i = 1:numberOfGenes
        current_gene_expression = human_expression(:,i);
        [expressionDist12, expressionDist13] = calcTripletExpresionDistance(current_gene_expression, tripletsIndices);
        result(i) = checkAgrementBetweenDistances2(expressionDist12, expressionDist13, distance12BiggerInTree);
        
%         random_ordering_of_expression = current_gene_expression(randomPerm);
%         random_result(i) = checkAgrementBetweenDistances2( random_human_expression(:,i) , human_gross_region_vec, treeDistances, tripletsIndices);
%         [randomExpressionDist12, randomExpressionDist13] = calcTripletExpresionDistance(random_ordering_of_expression, tripletsIndices);
        
        random_result(i,:) = checkAgrementBetweenDistances2(expressionDist12, expressionDist13, distance12BiggerInTree_random);
        
%         random_result(i,:) = checkAgrementBetweenDistances2_alt(current_gene_expression, randomPerm,tripletsIndices, distance12BiggerInTree);

        printPercentCounter(i,numberOfGenes);
%         fprintf('.');
    end

    
    
end

function [pairsIndices, treeDistances] = getTreeDistanceForAllPairs(region_vec, humanOntology, region_names)
    numberOfSamples = size(human_gross_region_vec,1);
    pairsIndices = nchoosek(1:numberOfSamples, 2);
    treeDistances = region_vec(pairsIndices);
end

function [expressionDist12, expressionDist13] = calcTripletExpresionDistance(expressionMatrix, tripletsIndices)
    [numberOfSamples, numberOfGenes] = size(expressionMatrix);
    expressionDist12 = abs(expressionMatrix(tripletsIndices(:,1),:) - expressionMatrix(tripletsIndices(:,2),:) );
    expressionDist13 = abs(expressionMatrix(tripletsIndices(:,1),:) - expressionMatrix(tripletsIndices(:,3),:) );
end

function [distance12, distance13] = calcTripletTreeDistance(distanceMatrix,region_vec, tripletsIndices)
    numberOfRegions = size(distanceMatrix,1);
    
    % find the index in the flat matrix.
    sampleRegionIndexInLongVec12 =  region_vec(tripletsIndices(:,1)) + (region_vec(tripletsIndices(:,2)) -1) * numberOfRegions;
    distance12 = distanceMatrix( sampleRegionIndexInLongVec12 );
    sampleRegionIndexInLongVec13 =  region_vec(tripletsIndices(:,1)) + (region_vec(tripletsIndices(:,3)) -1) * numberOfRegions;
    distance13 = distanceMatrix( sampleRegionIndexInLongVec13 );
end

function percent = checkAgrementBetweenDistances2(expressionDist12, expressionDist13, distance12BiggerInTree)

    [number_of_triplets, repeats] = size(distance12BiggerInTree);
    distance12BiggerInExpression = expressionDist12 > expressionDist13;
    S1 = sum(distance12BiggerInTree,1 );
    S2 = nan(1,repeats);
    for i = 1:repeats
        releventExp = distance12BiggerInExpression & distance12BiggerInTree(:,i) ;
        S2(i) = sum( releventExp ,1);
    end

%     releventExp = repmat(distance12BiggerInExpression,1,size(distance12BiggerInTree,2)) & distance12BiggerInTree ;
%     S1 = sum(distance12BiggerInTree,1 );
%     S2 = sum( releventExp ,1);

    percent = S2 ./ S1;

end

function percent = checkAgrementBetweenDistances2_alt(gene_expression, randomPerm, tripletsIndices, distance12BiggerInTree)


    random_ordering_of_expression = gene_expression(randomPerm);
    [expressionDist12_random, expressionDist13_random] = calcTripletExpresionDistance(random_ordering_of_expression, tripletsIndices);
    distance12BiggerInExpression_random = expressionDist12_random > expressionDist13_random;
    
    S1 = sum(distance12BiggerInTree,1 );
    S2 = sum( distance12BiggerInExpression_random(distance12BiggerInTree,:) ,1);
    percent = S2 ./ S1;

end


function result = checkAgrementBetweenDistances(expressionMatrix, region_vec, humanOntology, region_names,tripletsIndices)

    [numberOfSamples, numberOfGenes] = size(expressionMatrix);
    
     expressionDistances = abs(expressionMatrix(pairsIndices(:,1) ,: ) - expressionMatrix(pairsIndices(:,2) ,: ) );
%     expressionDistances = getExpressionDistanceAllPairs(expressionMatrix, pairsIndices);
%     [pairsIndices, treeDistances] = getTreeDistanceForAllPairs(human_gross_region_vec, humanOntology, region_names);
    
    n = numberOfSamples;
    distanceBiggerInTree = false( n*(n-1)*(n-2)/2 ,1);
    distanceBiggerInExpression = false( n*(n-1)*(n-2)/2 ,numberOfGenes);
%     indicesInCheck = zeros(n*(n-1)*(n-2)/2 ,3,'int8');
    j = 1;
    numberOfPairsWithoutOne = (n-1)*(n-2)/2;
    for i = 1:n
%         indicesWithout_i = [ 1: (i-1) , (i+1): numberOfSamples ];
%         pairsIndicesWithout_i = nchoosek(indicesWithout_i, 2);

        pairsIndicesWithout_i = pairsIndices( (pairsIndices(:,1) ~= i ) & (pairsIndices(:,2) ~= i) ,:);
        
        
        firstIndcies = [ones(size(pairsIndicesWithout_i,1),1)*i, pairsIndicesWithout_i(:,1)];
        firstIndcies = sort(firstIndcies,2);
        [~, firstDistance] = ismember( firstIndcies, pairsIndicesWithout_i,'rows');

        secondIndcies = [ones(size(pairsIndicesWithout_i,1),1)*i, pairsIndicesWithout_i(:,2)];
        secondIndcies = sort(secondIndcies,2);
        [~, secondDistance] = ismember( secondIndcies, pairsIndicesWithout_i,'rows');
        
        distanceBiggerInTreeTemp = treeDistances(firstDistance) > treeDistances(secondDistance);
        distanceBiggerInExpressionTemp = expressionDistances(firstDistance,:) > expressionDistances(secondDistance,:);
        
        distanceBiggerInTree(j: j + numberOfPairsWithoutOne-1) = distanceBiggerInTreeTemp;
        distanceBiggerInExpression(j: j + numberOfPairsWithoutOne-1,:) = distanceBiggerInExpressionTemp;
%         indicesInCheck(j:j + numberOfPairsWithoutOne-1,:) = [ones(size(pairsIndicesWithout_i,1),1)*i, pairsIndicesWithout_i];
        j = j + numberOfPairsWithoutOne;
    end

    S1 = sum(distanceBiggerInTree);
    S2 = sum( distanceBiggerInExpression(distanceBiggerInTree,:),1);
    
    result = (S1 ./ S2);
end

function randomIndicesWithNoRepeats = createRandomTriplets(n,d)
% create random d samples out of (n over 3) triplets

    randomIndices = randi(n,[d,3]);
    problamaticIndices12 = diff( randomIndices(:,[1,2]),1,2) == 0 ;
    problamaticIndices23 = diff( randomIndices(:,[2,3]),1,2) == 0 ;
    problamaticIndices13 = diff( randomIndices(:,[1,3]),1,2) == 0 ;


    randomIndicesWithNoRepeats =  randomIndices( ~problamaticIndices12 & ~problamaticIndices23 & ~problamaticIndices13 ,:);
    randomIndicesWithNoRepeats = unique(randomIndicesWithNoRepeats,'rows');
end

function [pairsIndices, expressionDistances] = getExpressionDistanceAllPairs(expressionMatrix,pairsIndices)
    [numberOfSamples, numberOfGenes] = size(expressionMatrix);
    pairsIndices = nchoosek(1:numberOfSamples, 2);
    expressionDistances = abs(expressionMatrix(pairsIndices(:,1) ,: ) - expressionMatrix(pairsIndices(:,2) ,: ) );
end