function treeMatrix = createBinaryTree(number_of_levels)

    numberOfNodes = 2^number_of_levels -1;
    treeMatrix = logical(sparse(numberOfNodes,numberOfNodes));
    
    [~, treeMatrix] = recursiveBuild(1, treeMatrix,1,number_of_levels);
end

function [lastOccuipedIndex, matrix] = recursiveBuild(index, matrix,level,max_level)
    
    lastOccuipedIndex = index;
    if (level < max_level)
        newNodeIndex = lastOccuipedIndex +1;
        matrix(index, newNodeIndex) = 1;
        [lastOccuipedIndex, matrix] = recursiveBuild(newNodeIndex, matrix,level +1,max_level);
        
        newNodeIndex = lastOccuipedIndex +1;
        matrix(index, newNodeIndex) = 1;
        [lastOccuipedIndex, matrix] = recursiveBuild(newNodeIndex, matrix,level +1,max_level);
    end
    
end