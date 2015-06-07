function zapala_ontology = build_zapala_ontlotgy()

    structureLabels = {'Grey matter';'Telencephalon';'Pallium';'Olfactory Bulbs';'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';'Hippocampus';'Hippocampus CA1';'Hippocampus CA3';'Dentate Gyrus';'Subpallium';'Amygdala';'Bed nucleus of the stria terminalis';'Striatum';'Diencephalon';'Thalamus';'Hypothalamus';'Mesencephalon';'Inferior Colliculi';'Superior Colliculi';'Metencephalon';'Cerebellum';'Pons';'Myelencephalon';'Medulla'};
    structureLabels = [cell(length(structureLabels),3) , structureLabels];
    dependencyMatrix = sparse(27,27);
    dependencyMatrix(1,2) = 1;
    dependencyMatrix(2,3) = 1;
    dependencyMatrix(3,4) = 1;
    dependencyMatrix(3,5) = 1;
    dependencyMatrix(5,6:8) = 1;
    dependencyMatrix(3,9) = 1;
    dependencyMatrix(9,10:12) = 1;
    dependencyMatrix(2,13) = 1;
    dependencyMatrix(13,14:16) = 1;
    dependencyMatrix(1,17 ) = 1;
    dependencyMatrix(17,18:19) = 1;
    dependencyMatrix(1,20) = 1;
    dependencyMatrix(20,21:22) = 1;
    dependencyMatrix(1,23) = 1;
    dependencyMatrix(23,24:25) = 1;
    dependencyMatrix(1,26) = 1;
    dependencyMatrix(26,27) = 1;
%     print_dependencies(dependencyMatrix,structureLabels(:,4));
    
    zapala_ontology.dependencyMatrix = dependencyMatrix;
    zapala_ontology.structureLabels = structureLabels;
    zapala_ontology.unDirectedDistanceMatrix = computeDistanceBetweenNodes(dependencyMatrix);
end

function print_dependencies(dependencyMatrix,list_of_regions)
    for i=2:size(dependencyMatrix,1)
       fprintf( '%s \t is child of \t %s\n' , list_of_regions{i} ,list_of_regions{find(dependencyMatrix(:,i)) });
    end
end

function unDirectedDistanceMatrix = computeDistanceBetweenNodes(dependencyMatrix)
    addpath('~/Projects/matlab_bgl')

    undirectedMatrix = dependencyMatrix + dependencyMatrix';
    directedDistanceMatrix = nan(size(dependencyMatrix));
    unDirectedDistanceMatrix = nan(size(dependencyMatrix));
    for i = 1:size(dependencyMatrix,1)
        [nodeDistance ~] = dijkstra_sp(dependencyMatrix,i);
        directedDistanceMatrix(:,i) = nodeDistance;
        
        [nodeDistance ~] = dijkstra_sp(undirectedMatrix,i);
        unDirectedDistanceMatrix(:,i) = nodeDistance;
    end
end