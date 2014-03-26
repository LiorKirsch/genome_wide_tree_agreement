function expresionVector = createTreeExpressionReverse(treeMatrix, region_vector)
    startNode = 1;
    startValue = 8;

    expresionVector = nan(size(region_vector,1),1);
    tree_expression = nan(size(treeMatrix,1),1);
    
    for i = 1:length(region_vector)
       current_region_index = region_vector(i); 
        
    end

    tree_expression(startNode) = startValue;
    [expresionVector, tree_expression] = recursiveExpression(expresionVector, treeMatrix, tree_expression, 1);
    
end

function [expresionVector, tree_expression] = recursiveExpression(expresionVector, treeMatrix, tree_expression, indexOfNode)

    inheretence_noise = 2;
    sample_noise = 0.5;
    
    
    if isnan(tree_expression(indexOfNode) )
        parentOfNode = find(treeMatrix(:,index));
        assert(length(parentOfNode) <= 1,'a node should have at most one parent');


        if isnan(tree_expression(parentOfNode) )
            [expresionVector, tree_expression] = recursiveExpression(expresionVector, treeMatrix, tree_expression, index);
        end

        parent_tree_expresion = tree_expression(parentOfNode);

        random_noise = rand( length(sonsOfNode) ,1)*2 -1;
        random_noise = random_noise * noise_amplitude;
        for i=1: length(sonsOfNode)
            sonIndex = sonsOfNode(i);
            sonValue = expresionVector(index) + random_noise(i);
            expresionVector( sonIndex ) = sonValue;
            expresionVector = recursiveExpression(expresionVector, treeMatrix, sonIndex);
        end
    end
    
    node_tree_expresion = tree_expression(indexOfNode);
    
end