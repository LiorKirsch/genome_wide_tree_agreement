function treeOrderingMainSimple(startIndex, finishIndex)

    addpath('/home/lab/gal/develop/matlab');
    addpath('~/Projects/general use functions/');

    
    addpath('/home/lab/lior/Projects/buildStructureOntology/');
    
    % Load the tree
    filename = '/home/lab/lior/Projects/buildStructureOntology/humanOntologyObject.mat';
    load(filename, 'humanOntology'); %#ok
    
    % Load the tree
    filename = '/home/lab/lior/Projects/buildStructureOntology/humanOntologyObjectV2.mat';
    humanOntologyV2 = load(filename, 'humanOntology'); %#ok
    humanOntologyV2 = humanOntologyV2.humanOntology;
    
    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);
%     fileName = 'results/simpleMeasurementGrossRegions';
% 
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6AllRegions',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
%     fileName = 'results/simpleMeasurementAllRegions';

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);
%     fileName = 'results/simpleMeasurementGrossCortexRegions';
    
    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',[]);
    fileName = 'results/simpleMeasurementAllCortexRegions';


    %normalize the data to have zero mean and unit variance
%     human_expression = human_expression - repmat(mean(human_expression,2),1,size(human_expression,2) );
%     human_expression = human_expression ./ repmat( std(human_expression,0,2),1,size(human_expression,2) );
%     developing_expression = developing_expression - repmat(mean(developing_expression,2),1,size(developing_expression,2) );
%     developing_expression = developing_expression ./ repmat( std(developing_expression,0,2),1,size(developing_expression,2) );
    %-----------------------------------------------------%

    medianGeneIndex = 942;   %not norm - COX7C
    medianGeneIndex = 18225; %Norm - LCE1E
    neurod1GeneIndex = find(strcmp('NEUROD1', human_gene_info.gene_symbols)); 

%     if any(strcmp(developing_gross_region_names, 'Hippocampus'))
%         developing_gross_region_names{ strcmp(developing_gross_region_names, 'Hippocampus') } = 'hippocampal formation';
%     end
    translated_developing_gross_region_names = translate_kang_region_names_to_ABA(developing_gross_region_names);
    
    [~,gross_region_indices_in_ontology] = ismember(human_gross_region_names, humanOntology.structureLabels(:,4) );
    [~,gross_developing_region_indices_in_ontology] = ismember(translated_developing_gross_region_names, humanOntologyV2.structureLabels(:,4) );
    
    human6TreeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology);
    developingTreeDistances = getTreeDistance(humanOntologyV2, gross_developing_region_indices_in_ontology);
    human6PhysicalDistances = squareform( pdist(physicalLocation,'euclidean') );
    
    samplesTreeDistance = getSamplesDistance(human_gross_region_vec, human6TreeDistances);
    samplesDevelopingTreeDistance = getSamplesDistance(developing_gross_region_vec, developingTreeDistances);
    
    
    identity_region_mapping = (1:length(human_gross_region_vec))';

%     drawTheValuesOfGene(neurod1GeneIndex, human_expression, samplesTreeDistance, human_gene_info);
%     drawTheValuesOfGene(medianGeneIndex, human_expression, samplesTreeDistance, human_gene_info);
    
    [human6Results,human6RandomResults] = compareExpressionToDistanceMetric(human_expression , samplesTreeDistance,startIndex, finishIndex);
    [brainspanResults,brainspanRandomResults] = compareExpressionToDistanceMetric(developing_expression, samplesDevelopingTreeDistance ,startIndex, finishIndex);

    human6PhysicalResults = nan(size(human6Results));
    human6PhysicalRandomResults = nan(size(human6RandomResults));
%     [human6PhysicalResults,human6PhysicalRandomResults] = compareExpressionToDistanceMetric(human_expression , human6PhysicalDistances,startIndex, finishIndex);
    

    fileName = sprintf('%s-%d-%d.mat',fileName, startIndex, finishIndex);
    indices = [startIndex, finishIndex];
    save(fileName, 'human6Results', 'human6RandomResults','human_gene_info','brainspanResults','brainspanRandomResults','developing_genes_info', 'human6PhysicalResults','human6PhysicalRandomResults','indices');
end

function [human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation] = drawXsamplesInRandom(numberOfSamples, human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation)
    [num_samples_human, num_genes] = size(human_expression);
    rand_indicies = randperm(num_samples_human);
    rand_indicies = rand_indicies(1:numberOfSamples);
    
    selected_samples = false(num_samples_human,1);
    selected_samples(rand_indicies) = true;
    
    human_expression = human_expression(selected_samples,:);
    human_samples2subjects = human_samples2subjects(selected_samples,:);
    physicalLocation = physicalLocation(selected_samples,:);
    human_gross_region_vec = human_gross_region_vec(selected_samples);    
end

function drawTheValuesOfGene(geneIndex, human_expression, distancesMatrix, human_gene_info)

    [numberOfSamples, numberOfGenes] = size(human_expression);

    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distancesMatrix(onlyUpperTri);

    current_gene_expression = human_expression(:,geneIndex);
    expression_distance_matrix = squareform( pdist(current_gene_expression,'euclidean') );
    onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
    
    spearmanScore = corr(onlyUpperExpressionMatrix, onlyUpperDistanceMatrix , 'type','Spearman');
    geneName = human_gene_info.gene_symbols{geneIndex};
    randomIndex = randi(length(onlyUpperDistanceMatrix), 10^5,1); % because there are too many dots to draw
    
    createFigure;
    hist2_alt( onlyUpperDistanceMatrix(randomIndex), onlyUpperExpressionMatrix(randomIndex),16,16);
    xlabel('Tree distance'); ylabel('Expression distance');
    fileName = sprintf('typical-%s-hist2.png', geneName);
    saveFigure(gcf, fileName, 'png');
    
    createFigure;
    onlyUpperDistanceMatrixWithNoise = onlyUpperDistanceMatrix + (rand(size(onlyUpperDistanceMatrix))*0.2 - 0.1 );
    scatter(onlyUpperDistanceMatrixWithNoise(randomIndex), onlyUpperExpressionMatrix(randomIndex),5,'filled' );
    title( sprintf('%s   spearman: %g' , geneName, spearmanScore ));
    xlabel('Tree distance','Fontsize',20);     ylabel('Expression distance','Fontsize',20);
    fileName = sprintf('typical-%s.png', geneName);
    saveFigure(gcf, fileName, 'png');
    
    createFigure;
    transparentScatter(  onlyUpperDistanceMatrixWithNoise(randomIndex), onlyUpperExpressionMatrix(randomIndex) ,0.01,0.3);
    xlabel('Tree distance','Fontsize',20);     ylabel('Expression distance','Fontsize',20);
    fileName = sprintf('typicalTrans-%s.png', geneName);
    saveFigure(gcf, fileName, 'png');

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
    
    repeat_random = 100;
    randomPerm = nan(numberOfSamples , repeat_random);
    for i = 1:repeat_random
        randomPerm(:,i) = randperm(numberOfSamples);
    end
    random_human_expression = human_expression(randperm(numberOfSamples),:);
    
    

    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distancesMatrix(onlyUpperTri);
    
    onlyUpperRandomDistanceMatrix =  nan(length(onlyUpperDistanceMatrix) ,repeat_random);
    for i = 1:repeat_random
        randomDistanceMatrix = distancesMatrix( randomPerm(:,i), randomPerm(:,i) );
        onlyUpperRandomDistanceMatrix(:,i) = randomDistanceMatrix(onlyUpperTri);
    end
    
    tiedrank_randomTreeDistance = tiedrank( onlyUpperRandomDistanceMatrix);
    tiedrank_treeDistance = tiedrank( onlyUpperDistanceMatrix);
    
  
    [result, random_result] = calcCorrInForAllGenes(human_expression, tiedrank_treeDistance, tiedrank_randomTreeDistance, startIndex, finishIndex);

    
    
end

function [result, random_result] = calcCorrInForAllGenes(human_expression, tiedrank_treeDistance, tiedrank_randomTreeDistance, startIndex, finishIndex)
    
    [numberOfSamples, numberOfGenes] = size(human_expression);
    repeat_random = size(tiedrank_randomTreeDistance,2);
    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    numberPerSlice = 25;
    
    result = nan(numberOfGenes,1);
    random_result = nan(numberOfGenes,repeat_random);
    
 
    lastIndex = min(finishIndex, numberOfGenes);
    fprintf('finished preprocessing, entering gene loop ----- ');
    fprintf('progress:    ');
    
    indicesForRandom = round(linspace(1,repeat_random,round(repeat_random/numberPerSlice) ));

    for j = 1:(length(indicesForRandom) -1)
        curr_random_result = nan(numberOfGenes,(1+ indicesForRandom(j+1) - indicesForRandom(j)) ) ;
        sliced_tiedrank_randomTreeDistance = tiedrank_randomTreeDistance(:, indicesForRandom(j): indicesForRandom(j+1));
        
        parfor i = startIndex: lastIndex
            current_gene_expression = human_expression(:,i);
            expression_distance_matrix = squareform( pdist(current_gene_expression,'euclidean') );
            onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
            tiedrank_expression = tiedrank(onlyUpperExpressionMatrix);
%             tiedrank_expression = single( tiedrank_expression) ;
            result(i) = corr(tiedrank_expression, tiedrank_treeDistance , 'type','Pearson');
            curr_random_result(i,:) = corr(tiedrank_expression, sliced_tiedrank_randomTreeDistance , 'type','Pearson');

        end
        random_result(:, indicesForRandom(j): indicesForRandom(j+1))  = curr_random_result;
        fprintf( '\b\b\b%2d%%', round(j/(length(indicesForRandom) -1)*100) );
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

