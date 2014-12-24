function treeOrderingMainSimpleSingleSubject(startIndex, finishIndex)

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
    

%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);
%     fileName = 'results/simpleMeasurementGrossRegions';

    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6AllRegions',[]);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
    fileName = 'results/simpleMeasurementAllRegions';

    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangAllRegions',[]);
%     fileName = 'results/simpleMeasurement16Regions';
    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6Cortex',[]);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortex',[]);
%     fileName = 'results/simpleMeasurementGrossCortexRegions';
    
%     [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6CortexAllRegions', []);
%     [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kangCortexAllRegions',[]);
%     fileName = 'results/simpleMeasurementAllCortexRegions';



     [human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation] = drawXsamplesInRandom( 500 , human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation);
    
    %normalize the data to have zero mean and unit variance
%     human_expression = human_expression - repmat(mean(human_expression,2),1,size(human_expression,2) );
%     human_expression = human_expression ./ repmat( std(human_expression,0,2),1,size(human_expression,2) );
    %-----------------------------------------------------%

    translated_developing_gross_region_names = translate_kang_region_names_to_ABA(developing_gross_region_names);

    [~,gross_region_indices_in_ontology] = ismember(human_gross_region_names, humanOntology.structureLabels(:,4) );
    [~,gross_developing_region_indices_in_ontology] = ismember(translated_developing_gross_region_names, humanOntologyV2.structureLabels(:,4) );

    
    human6TreeDistances = getTreeDistance(humanOntology, gross_region_indices_in_ontology);
    developingTreeDistances = getTreeDistance(humanOntologyV2, gross_developing_region_indices_in_ontology);

    samplesTreeDistance = getSamplesDistance(human_gross_region_vec, human6TreeDistances);
    samplesDevelopingTreeDistance = getSamplesDistance(developing_gross_region_vec, developingTreeDistances);

    
    identity_region_mapping = (1:length(human_gross_region_vec))';
 
    num_subjects_dev = size(developing_samples2subjects,2);
    brainspanResults = cell(num_subjects_dev,1);
    brainspanRandomResults = cell(num_subjects_dev,1);
    for i = 1:num_subjects_dev
        current_subject_samples = developing_samples2subjects(:,i) ;
        [brainspanResults{i},brainspanRandomResults{i}] = compareExpressionToDistanceMetric(developing_expression(current_subject_samples,:), samplesDevelopingTreeDistance(current_subject_samples,current_subject_samples) ,startIndex, finishIndex);
        fprintf('done developing subject %d\n',i);
    end
    
    num_subjects = size(human_samples2subjects,2);
    human6Results = cell(num_subjects,1);
    human6RandomResults = cell(num_subjects,1);
    for i = 1:num_subjects
        current_subject_samples = human_samples2subjects(:,i) ;
%         figure('name',sprintf('distance distribution subject %d', i));
%         drawDistancesDistribution( samplesTreeDistance(current_subject_samples,current_subject_samples) );
        [human6Results{i},human6RandomResults{i}] = compareExpressionToDistanceMetric(human_expression(current_subject_samples,:) , samplesTreeDistance(current_subject_samples,current_subject_samples),startIndex, finishIndex);
        fprintf('done subject %d\n',i);
    end
  
  
     

    fileName = sprintf('%s-%d-%d-subjects.mat',fileName, startIndex, finishIndex);
    indices = [startIndex, finishIndex];
    save(fileName, 'human6Results', 'human6RandomResults','brainspanResults','brainspanRandomResults','human_gene_info', 'developing_genes_info','indices','-v7.3');
end

function drawDistancesDistribution(distanceMatrix)
    allDistances = distanceMatrix(:);
    bins = 0:1:(max(allDistances) +1)  ;
    countInBins = histc(allDistances,bins -0.5) / length( allDistances);
    bar(bins, countInBins);
    xlabel('distance');
    ylabel('normalized count');
end

function [human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation] = drawXsamplesInRandom(numberOfSamples, human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation)
    [num_samples_human, num_genes] = size(human_expression);
%     s = RandStream('mt19937ar','Seed','shuffle');
%     RandStream.setGlobalStream(s);
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
    
    repeat_random = 5000;
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
    fprintf('\n');
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

