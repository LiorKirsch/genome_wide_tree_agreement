function drawGeneAcrossRegions(gene_sym)
    
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
    

    [human_expression, human_gross_region_vec, human_gene_info, ~,...
        human_gross_region_names, ~] = load_expression_and_regions('human6GrossRegions',[]);
  
    region_color = nan(length(human_gross_region_names),3);
    for i =1 :length(human_gross_region_names)
        region_color(i,:) = humanOntology.getColorByRegionName(human_gross_region_names{i});
    end
    
    
    
    for i =1:length(gene_sym)
        fprintf('drawing expression bars for gene %s\n', gene_sym{i});
        geneIndex = find(strcmp(gene_sym{i}, human_gene_info.gene_symbols));


        expression_of_gene = human_expression(:,geneIndex);
        mean_exp = accumarray(human_gross_region_vec, expression_of_gene, [], @mean);
        std_exp = accumarray(human_gross_region_vec, expression_of_gene, [], @std);

        figure;
        drawBars(mean_exp,std_exp,human_gross_region_names,region_color,'mean expression');

        title(gene_sym{i});
        fileName = sprintf('figures/bars_brain_%s', gene_sym{i});
        saveFigure(gcf, [fileName, '.png'], 'png');
%         saveFigure(gcf, [fileName, '.tiff'], 'tiff');
        saveFigure(gcf, [fileName, '.eps'], 'eps');
    end
end