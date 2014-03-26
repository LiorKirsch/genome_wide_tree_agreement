function drawExpressionOfGene()

    filename = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability','humanOntologyObject.mat');
    load(filename, 'humanOntology'); %#ok

    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names, physicalLocation] = load_expression_and_regions('human6',[]);
    [developing_expression, developing_gross_region_vec, developing_genes_info, developing_samples2subjects, developing_gross_region_names] = load_expression_and_regions('kang',[]);

     if any(strcmp(developing_gross_region_names, 'Hippocampus'))
        developing_gross_region_names{ strcmp(developing_gross_region_names, 'Hippocampus') } = 'hippocampal formation';
    end
    
    human_gross_region_colors = humanOntology.getColorByRegionName(human_gross_region_names);
    dev_human_gross_region_colors = humanOntology.getColorByRegionName(developing_gross_region_names);

    gene_name = 'SYNPO';
    gene_index = strcmp(human_gene_info.gene_symbols , gene_name);
    gene_index_dev = strcmp(developing_genes_info.gene_symbols , gene_name);
    
    
    figure('name',sprintf('Human6 expression - %s', gene_name) );
    gene_expression = human_expression(:,gene_index);
    meanExpression = accumarray(human_gross_region_vec,gene_expression,[],@(x)mean(x,1));
    stdExpression = accumarray(human_gross_region_vec,gene_expression,[],@(x)std(x,1,1));
    
%     plot(1:length(meanExpression), meanExpression);
%     xticklabel_rotate(1:length(meanExpression),45,human_gross_region_names);
    drawBars(meanExpression,stdExpression,human_gross_region_names,human_gross_region_colors,'mean expression');
    
    
    figure('name',sprintf('Developing human expression - %s', gene_name) );
    gene_expression = developing_expression(:,gene_index_dev);
    meanExpression = accumarray(developing_gross_region_vec,gene_expression,[],@(x)mean(x,1));
    stdExpression = accumarray(developing_gross_region_vec,gene_expression,[],@(x)std(x,1,1));
    drawBars(meanExpression,stdExpression,developing_gross_region_names,dev_human_gross_region_colors,'mean expression');
end