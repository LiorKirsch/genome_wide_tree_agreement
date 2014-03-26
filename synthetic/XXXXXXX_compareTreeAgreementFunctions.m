function compareTreeAgreementFunctions()

    conf = configureation();
    
    %============== Load the tree ==================
    addpath('~/Projects/individual variability');
    load('~/Projects/individual variability/humanOntologyObject.mat');
    treeMatrix = humanOntology.dependencyMatrix;
    
%     %======= create a synthetic binary tree ========
%     treeMatrix = createBinaryTree(10);
    
    samples_region_vector = 1:length(treeMatrix) ; % one sample per area
    expressionVector = createTreeExpression(treeMatrix, conf);

    unDirectedDistanceMatrix = computeDistanceBetweenNodes(double(treeMatrix));
    
    result = agreementUsingCorr(expressionVector, unDirectedDistanceMatrix);
    
%     plotResultsForDimentions(treeMatrix, conf)

end


function result = agreementUsingCorr(expressionVector, distanceMatrix)
    numberOfSamples = size(expressionVector,1);
    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distanceMatrix(onlyUpperTri);
   
    expression_distance_matrix = squareform( pdist(expressionVector,'euclidean') );
    onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
    result = corr(onlyUpperExpressionMatrix, onlyUpperDistanceMatrix , 'type','Spearman');
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

function plotResultsForDimentions(treeMatrix, conf)
    unDirectedDistanceMatrix = computeDistanceBetweenNodes(double(treeMatrix));

    dim = 2.^(0:10);
    results = nan(size(dim));
    for i= 1:length(dim)
        conf.expressionDimention = dim(i);
        expressionVector = createTreeExpression(treeMatrix, conf);
        results(i) = agreementUsingCorr(expressionVector, unDirectedDistanceMatrix);
    end
    semilogx(dim,results);
end