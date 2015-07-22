function drawTheValuesOfGeneMain()

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

%     [human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation] = drawXsamplesInRandom(length(developing_gross_region_vec), human_expression, human_gross_region_vec, human_samples2subjects, physicalLocation);
    
    %normalize the data to have zero mean and unit variance
%     human_expression = human_expression - repmat(mean(human_expression,2),1,size(human_expression,2) );
%     human_expression = human_expression ./ repmat( std(human_expression,0,2),1,size(human_expression,2) );
%     developing_expression = developing_expression - repmat(mean(developing_expression,2),1,size(developing_expression,2) );
%     developing_expression = developing_expression ./ repmat( std(developing_expression,0,2),1,size(developing_expression,2) );
    %-----------------------------------------------------%

    medianGeneIndex = 942;   %not norm - COX7C
    medianGeneIndex = 18225; %Norm - LCE1E
    neurod1GeneIndex = find(strcmp('NEUROD1', human_gene_info.gene_symbols)); 
    synpoGeneIndex = find(strcmp('SYNPO', human_gene_info.gene_symbols)); 
    

    translated_developing_gross_region_names = translate_kang_region_names_to_ABA(developing_gross_region_names);
    
    [~,gross_region_indices_in_ontology] = ismember(human_gross_region_names, humanOntology.structureLabels(:,4) );

    human6TreeDistances = humanOntology.unDirectedDistanceMatrix(gross_region_indices_in_ontology,gross_region_indices_in_ontology);
    samplesTreeDistance = getSamplesDistance(human_gross_region_vec, human6TreeDistances);
    
    drawTheValuesOfGene(synpoGeneIndex, human_expression, samplesTreeDistance, human_gene_info);
    drawTheValuesOfGene(neurod1GeneIndex, human_expression, samplesTreeDistance, human_gene_info);
%     drawTheValuesOfGene(medianGeneIndex, human_expression, samplesTreeDistance, human_gene_info);
  
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
    
    maxTreeDistance = max(onlyUpperDistanceMatrix(randomIndex));
    hist2_alt( onlyUpperDistanceMatrix(randomIndex), onlyUpperExpressionMatrix(randomIndex),maxTreeDistance,16);
    set(gca,'FontSize',20);
    xlabel('Ontology distance'); ylabel('Expression distance');
    fileName = sprintf('figures/2dhist-%s-corr-%.2g', geneName, spearmanScore);
    saveFigure(gcf, [fileName, '.png'], 'png');
    saveFigure(gcf, [fileName, '.tiff'], 'tiff');
    saveFigure(gcf, [fileName, '.eps'], 'eps');

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