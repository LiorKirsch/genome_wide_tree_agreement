function treeOrderingMainSimple(startIndex, finishIndex)

    addpath('/home/lab/gal/develop/matlab');
    addpath('~/Projects/general use functions/');
    addpath('/home/lab/lior/Projects/buildStructureOntology/');
   
    kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};


% dataset_name = 'human6AllRegions'; % ======= human6 - all regions =======
% dataset_name = 'human6GrossRegions'; % ======= human6 - gross regions =======
% dataset_name = 'human6Cortex'; % ======= human6 - cortex gross regions=======
% dataset_name = 'human6CortexAllRegions'; % ======= human6 - cortex all regions=======

% dataset_name = 'kangAllRegions'; % ======= Kang - all regions =======
% dataset_name = 'kangCortexAllRegions'; % ======= Kang - cortex all regions =======
% dataset_name = 'kangGrossRegions'; % ======= Kang - gross regions =======
% dataset_name = 'kangCortex'; % ======= Kang - cortex gross regions =======

dataset_name = 'zapalaMouse'; % ======= Zapala mouse =======


fileName = fullfile('results',dataset_name);
[expression, gross_region_vec, gene_info, samples2subjects, gross_region_names, physicalLocation] = load_expression_and_regions(dataset_name,[]);


% % ------- homology groups insead of genes -------
% addpath('/cortex/code/matlab/homologous_gene_mapping/');
% [~, ~, human_gene_info, ~, ~, ~] = load_expression_and_regions('human6AllRegions',[]);
% [~, ~, mouse_gene_info, ~, ~, ~] = load_expression_and_regions('zapalaMouse',[]);
% [gene_to_group_matrix_mouse, gene_to_group_matrix_human, homologous_group_id] = gene_to_homolog_group('mouse_laboratory','human', mouse_gene_info.gene_symbols, 'symbol',human_gene_info.gene_symbols,'symbol',false);
% switch dataset_name
%     case {'human6AllRegions','human6GrossRegions','human6Cortex','human6CortexAllRegions','kangAllRegions','kangCortexAllRegions','kangGrossRegions','kangCortex'}
%         expression = get_group_expression(gene_to_group_matrix_human, expression);
%     case 'zapalaMouse'
%         expression = get_group_expression(gene_to_group_matrix_mouse, expression);
% end
% gene_info.gene_symbols = homologous_group_id;
% gene_info.entrez_ids = homologous_group_id;
% fileName = [fileName,'_homologs'];
% expression = single(expression);


% -------normalize the data to have zero mean and unit variance-------
%     expression = expression - repmat(mean(expression,2),1,size(expression,2) );
%     expression = expression ./ repmat( std(expression,0,2),1,size(expression,2) );
    %-----------------------------------------------------%


    % ================== load ontology ==================
    switch dataset_name
        case {'human6AllRegions','human6GrossRegions','human6Cortex','human6CortexAllRegions'}
            selected_ontology = load('/home/lab/lior/Projects/buildStructureOntology/humanOntologyObject.mat', 'humanOntology'); 
            selected_ontology = selected_ontology.humanOntology;
            
        case {'kangAllRegions','kangCortexAllRegions','kangGrossRegions','kangCortex'}
            selected_ontology = load('/home/lab/lior/Projects/buildStructureOntology/humanOntologyObjectV2.mat', 'humanOntology'); 
            selected_ontology = selected_ontology.humanOntology;

        case 'zapalaMouse'
            selected_ontology = build_zapala_ontlotgy();
    end


    % ============= translate kang regions to ABA ========
    switch dataset_name
        case {'kangAllRegions','kangCortexAllRegions'}
            gross_region_names = translate_kang_region_names_to_ABA(gross_region_names);
        case {'kangGrossRegions','kangCortex'}
            gross_region_names{ strcmp(gross_region_names, 'Hippocampus') } = 'hippocampal formation';
    end

    
    [~,gross_region_indices_in_ontology] = ismember(gross_region_names, selected_ontology.structureLabels(:,4) );
    
    treeDistances = getTreeDistance(selected_ontology, gross_region_indices_in_ontology);
    samplesTreeDistance = getSamplesDistance(gross_region_vec, treeDistances);

%     selecedGeneSymbol = 'NEUROD1'; %'FEZF2','COX7C','LCE1E'
%     selectedGeneIndex = find(strcmp(selecedGeneSymbol, gene_info.gene_symbols)); 
%     drawTheValuesOfGene(selectedGeneIndex, expression, samplesTreeDistance, gene_info);

   
%     fprintf('BRO for average expression of genes, ');
%     [BRO_score_for_gene_average,~] = compareExpressionToDistanceMetric(mean(expression,2) , samplesTreeDistance,startIndex, finishIndex);
%     fprintf('BRO for %s, ',selecedGeneSymbol);
%     [BRO_score_for_selected_gene,random_BRO_score_for_selected_gene] = compareExpressionToDistanceMetric( expression(:,selectedGeneIndex) , samplesTreeDistance,startIndex, finishIndex);
%     figure('Name',sprintf('%s perm dist',selecedGeneSymbol)); hold on; hist(random_BRO_score_for_selected_gene,100) ;stem(BRO_score_for_selected_gene, 10,'r'); legend(sprintf('%s permuted',selecedGeneSymbol),selecedGeneSymbol);  hold off;
    
    fprintf('computing BRO score %s, ', dataset_name);
    [results,randomResults] = compareExpressionToDistanceMetric(expression , samplesTreeDistance,startIndex, finishIndex);
   

    fileName = sprintf('%s-%d-%d.mat',fileName, startIndex, finishIndex);
    indices = [startIndex, finishIndex];
    save(fileName, 'results', 'randomResults','gene_info','indices','-v7.3');
end

function group_expression = get_group_expression(gene_to_group_matrix, gene_expression)

    num_of_genes_in_group = sum(gene_to_group_matrix,1);
    assert( all(num_of_genes_in_group > 0), 'every group should have atleast one gene');
    
    group_expression = gene_expression * gene_to_group_matrix;
    group_expression = group_expression * diag( 1./ num_of_genes_in_group );
    
end

function [expression, gross_region_vec, samples2subjects, physicalLocation] = drawXsamplesInRandom(numberOfSamples, expression, gross_region_vec, samples2subjects, physicalLocation)
    [num_samples_human, num_genes] = size(expression);
    s = RandStream('mt19937ar','Seed','shuffle');
    RandStream.setGlobalStream(s);
    rand_indicies = randperm(num_samples_human);
    rand_indicies = rand_indicies(1:numberOfSamples);
    
    selected_samples = false(num_samples_human,1);
    selected_samples(rand_indicies) = true;
    
    expression = expression(selected_samples,:);
    samples2subjects = samples2subjects(selected_samples,:);
    physicalLocation = physicalLocation(selected_samples,:);
    gross_region_vec = gross_region_vec(selected_samples);    
end

function drawTheValuesOfGene(geneIndex, expression, distancesMatrix, gene_info)

    [numberOfSamples, numberOfGenes] = size(expression);

    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distancesMatrix(onlyUpperTri);

    current_gene_expression = expression(:,geneIndex);
    expression_distance_matrix = squareform( pdist(current_gene_expression,'euclidean') );
    onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
    
    spearmanScore = corr(onlyUpperExpressionMatrix, onlyUpperDistanceMatrix , 'type','Spearman');
    geneName = gene_info.gene_symbols{geneIndex};
    randomIndex = randi(length(onlyUpperDistanceMatrix), 10^5,1); % because there are too many dots to draw
    
    createFigure;
    hist2_alt( onlyUpperDistanceMatrix(randomIndex), onlyUpperExpressionMatrix(randomIndex),nan,16);
    xlabel('Tree distance'); ylabel('Expression distance');
    fileName = sprintf('figures/typical-%s-hist2', geneName);
    saveFigure(gcf, fileName, 'png');
    saveFigure(gcf, fileName, 'eps');
    
%     createFigure;
%     onlyUpperDistanceMatrixWithNoise = onlyUpperDistanceMatrix + (rand(size(onlyUpperDistanceMatrix))*0.2 - 0.1 );
%     scatter(onlyUpperDistanceMatrixWithNoise(randomIndex), onlyUpperExpressionMatrix(randomIndex),5,'filled' );
%     title( sprintf('%s   spearman: %g' , geneName, spearmanScore ));
%     xlabel('Tree distance','Fontsize',20);     ylabel('Expression distance','Fontsize',20);
%     fileName = sprintf('typical-%s.png', geneName);
%     saveFigure(gcf, fileName, 'png');
%     
%     createFigure;
%     transparentScatter(  onlyUpperDistanceMatrixWithNoise(randomIndex), onlyUpperExpressionMatrix(randomIndex) ,0.01,0.3);
%     xlabel('Tree distance','Fontsize',20);     ylabel('Expression distance','Fontsize',20);
%     fileName = sprintf('typicalTrans-%s.png', geneName);
%     saveFigure(gcf, fileName, 'png');

end

function treeDistances = getTreeDistance(selected_ontology, gross_region_indices_in_ontology)
%     treeDistances = selected_ontology.meanDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
    treeDistances = selected_ontology.unDirectedDistanceMatrix(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = selected_ontology.longDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
%   treeDistances = selected_ontology.shortDistanceToParent(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
end

function [result, random_result] = compareExpressionToDistanceMetric( expression, distancesMatrix ,startIndex, finishIndex)
    

    
    [numberOfSamples, numberOfGenes] = size(expression);
    
    rand('state', 1221);
    
    repeat_random = 100;
    numberPerBatch = repeat_random;
    randomPerm = nan(numberOfSamples , repeat_random);
    for i = 1:repeat_random
        randomPerm(:,i) = randperm(numberOfSamples);
    end
    random_human_expression = expression(randperm(numberOfSamples),:);
    
    
    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
    onlyUpperDistanceMatrix = distancesMatrix(onlyUpperTri);
    tiedrank_treeDistance = tiedrank( onlyUpperDistanceMatrix);    
    
    tiedrank_treeFullDistance = zeros( size(distancesMatrix) );
    tiedrank_treeFullDistance(onlyUpperTri) = tiedrank_treeDistance;  
    tiedrank_treeFullDistance = tiedrank_treeFullDistance + tiedrank_treeFullDistance';
    
    tiedrank_treeFullDistance = single(tiedrank_treeFullDistance);
    
    result = calcCorrInForAllGenes(expression, tiedrank_treeDistance, startIndex, finishIndex);
    
    random_result = nan(numberOfGenes,repeat_random);
    
 
    fprintf('finished preprocessing, entering gene loop ----- ');
    fprintf('progress:    ');
    
    indicesForRandom = round(linspace(1,repeat_random,round(repeat_random/numberPerBatch) +1 ));

    for j = 1:(length(indicesForRandom) -1)
        batchFirstIndex = indicesForRandom(j);
        batchLastIndex = indicesForRandom(j+1);
        batchSize = (1+ batchLastIndex - batchFirstIndex);
        onlyUpperRandomDistanceMatrix =  nan(length(onlyUpperDistanceMatrix) ,batchSize);
        tiedrank_randomTreeDistance =  nan(length(onlyUpperDistanceMatrix) ,batchSize);
        
%         for m = 1:batchSize
%             currentRandIndex = m + batchFirstIndex -1;
%             currentPerm = randomPerm(:,currentRandIndex);
%             randomDistanceMatrix = distancesMatrix( currentPerm, currentPerm );
%             onlyUpperRandomDistanceMatrix(:,m) = randomDistanceMatrix(onlyUpperTri);
%         end
%         tiedrank_randomTreeDistance = tiedrank( onlyUpperRandomDistanceMatrix);
        
        for m = 1:batchSize
            currentRandIndex = m + batchFirstIndex -1;
            currentPerm = randomPerm(:,currentRandIndex);
            permTiedRankDistances = tiedrank_treeFullDistance( currentPerm, currentPerm );
            tiedrank_randomTreeDistance(:,m) = permTiedRankDistances(onlyUpperTri);
        end

        curr_random_result = calcCorrInForAllGenes(expression, tiedrank_randomTreeDistance, startIndex, finishIndex);
    
        random_result(:, indicesForRandom(j): indicesForRandom(j+1))  = curr_random_result;
        fprintf( '\b\b\b%2d%%', round(j/(length(indicesForRandom) -1)*100) );
        
    end 
    fprintf('\n');
    
end

function result = calcCorrInForAllGenes(expression, tiedrank_treeDistance, startIndex, finishIndex)
    
    [numberOfSamples, numberOfGenes] = size(expression);
    repeats = size(tiedrank_treeDistance,2);
    onlyUpperTri = triu(true(numberOfSamples,numberOfSamples),1 );
     
    result = nan(numberOfGenes,repeats);
    
 
    lastIndex = min(finishIndex, numberOfGenes);
    
    parfor i = startIndex: lastIndex
        current_gene_expression = expression(:,i);
        expression_distance_matrix = squareform( pdist(current_gene_expression,'euclidean') );
        onlyUpperExpressionMatrix = expression_distance_matrix(onlyUpperTri);
        tiedrank_expression = tiedrank(onlyUpperExpressionMatrix);
%             tiedrank_expression = single( tiedrank_expression) ;
        result(i,:) = corr(tiedrank_expression, tiedrank_treeDistance , 'type','Pearson');
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

