function compareTreeAgreementFunctions()

    conf = configureation();
    
    %============== Load the tree ==================
    addpath('~/Projects/individual variability');
    addpath('~/Projects/general use functions/');
    addpath('~/Projects/genome_wide_agreement/');
    load('~/Projects/individual variability/humanOntologyObject.mat');
    treeMatrix = humanOntology.dependencyMatrix;
    [~, human_gross_region_vec, ~, ~, ~, ~] = load_expression_and_regions('human6AllRegions',[]);

    
%     %======= create a synthetic binary tree ========
%     treeMatrix = createBinaryTree(10);
    
    tree_node_expression = createTreeExpression(treeMatrix, conf);
    unDirectedDistanceMatrix = computeDistanceBetweenNodes(double(treeMatrix));
    resultForRegions = agreementUsingCorr(tree_node_expression, unDirectedDistanceMatrix);    
%     plotResultsForDimentions(treeMatrix, (1:size(treeMatrix,1))', conf)

    sample_expression  = create_region_expression( tree_node_expression, human_gross_region_vec, conf);
    distance_for_samples = distanceForSamplesUsingDistanceMatrix(unDirectedDistanceMatrix, human_gross_region_vec);
   
    resultForSamples = agreementUsingCorr(sample_expression, distance_for_samples);
    
%     plotResultsForDimentions(treeMatrix, human_gross_region_vec, conf);
end

function output = distanceForSamplesUsingDistanceMatrix(distanceMatrix, sample_region_index)
    number_of_samples = size(sample_region_index,1);
    output = nan(number_of_samples);
    
    for i = 1:number_of_samples
        sample_i_region = sample_region_index(i);
        output(i,:) = distanceMatrix(sample_i_region, sample_region_index);
    end
end


function result = agreementUsingCorr(expressionVector, tree_distances_sample_matrix)
% computes the correlation of the tree distance and the expression distance 
% for all pairs of samples

    numberOfSamples = size(expressionVector,1);
    assert( numberOfSamples == size(tree_distances_sample_matrix,1) && numberOfSamples == size(tree_distances_sample_matrix,2) ,'tree_distances_sample_matrix should be a n x n matrix');
    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = tree_distances_sample_matrix(onlyUpperTri);
   
    expression_distance_matrix = squareform( pdist(expressionVector,'euclidean') );
    onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
    result = corr(onlyUpperExpressionMatrix, onlyUpperDistanceMatrix , 'type','Spearman');
end


function result = agreementUsingTriplets(expressionVector, distanceMatrix)
% draw triplets from the data and computes the agreement

    
end

function directedDistences = getDistances(dependencyMatrix)
    X = dependencyMatrix;
    directedDistences = ( inv(eye(size(X)) - X) )^2*X;
    directedDistences(directedDistences ==0) = inf;
    directedDistences(logical(eye(size(X)))) = 0;
end

function unDirectedDistanceMatrix = computeDistanceBetweenNodes(dependecyMatrix)
    addpath('~/Projects/matlab_bgl')
  
    undirectedMatrix = dependecyMatrix + dependecyMatrix';
    directedDistanceMatrix = nan(size(dependecyMatrix));
    unDirectedDistanceMatrix = nan(size(dependecyMatrix));
    for i = 1:size(dependecyMatrix,1)
        [nodeDistance, ~] = dijkstra_sp(dependecyMatrix,i);
        directedDistanceMatrix(:,i) = nodeDistance;
        
        [nodeDistance, ~] = dijkstra_sp(undirectedMatrix,i);
        unDirectedDistanceMatrix(:,i) = nodeDistance;
        
    end
end

function plotResultsForDimentions(treeMatrix, sample_region_index, conf)
    unDirectedDistanceMatrix = computeDistanceBetweenNodes(double(treeMatrix));
    tree_distances = distanceForSamplesUsingDistanceMatrix(unDirectedDistanceMatrix, sample_region_index);

    dim = 2.^(0:10);
    dim = [1,2];
    repeat = 100;
    results = nan(length(dim),repeat);

    for i= 1:length(dim)
        conf.expressionDimention = dim(i);
        for j = 1:repeat
            tree_node_expression = createTreeExpression(treeMatrix, conf);
            sample_expression  = create_region_expression( tree_node_expression, sample_region_index, conf);

            results(i,j) = agreementUsingCorr(sample_expression, tree_distances);
        end
    end
    createFigure();
    h = plotp(dim,results');
    set(gca,'xscale','log');
    xlabel('dimention');
    ylabel('agreement score');
end