function expresionVector = createTreeExpression(treeMatrix, conf)

    expresionVector = zeros(size(treeMatrix,1),conf.expressionDimention);
    expresionVector( conf.root_node_index ,:) = conf.root_node_value ;
    expresionVector = forwardExpression(expresionVector, treeMatrix, 1, conf);
    
end

function expresionVector = forwardExpression(expresionVector, treeMatrix, index, conf)

    sonsOfNode = find(treeMatrix(index,:));
    
    random_noise = rand( length(sonsOfNode) ,conf.expressionDimention)*2 -1;
    random_noise = random_noise * conf.inheretence_noise;
    for i=1: length(sonsOfNode)
        sonIndex = sonsOfNode(i);
        sonValue = expresionVector(index,:) + random_noise(i,:);
        expresionVector( sonIndex ,:) = sonValue;
        expresionVector = forwardExpression(expresionVector, treeMatrix, sonIndex, conf);
    end
end