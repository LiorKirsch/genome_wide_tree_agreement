addpath('~/Projects/general use functions/');
addpath('/home/lab/lior/Projects/buildStructureOntology/');
   

dataset_name = 'human6AllRegions'; % ======= human6 - all regions =======
% dataset_name = 'human6GrossRegions'; % ======= human6 - gross regions =======
% dataset_name = 'human6Cortex'; % ======= human6 - cortex gross regions=======
% dataset_name = 'human6CortexAllRegions'; % ======= human6 - cortex all regions=======

[human_expression, ~, human_gene_info, ~, ~, ~] = load_expression_and_regions(dataset_name,[]);

% dataset_name = 'kangAllRegions'; % ======= Kang - all regions =======
% dataset_name = 'kangCortexAllRegions'; % ======= Kang - cortex all regions =======
% dataset_name = 'kangGrossRegions'; % ======= Kang - gross regions =======
% dataset_name = 'kangCortex'; % ======= Kang - cortex gross regions =======

dataset_name = 'zapalaMouse'; % ======= Zapala mouse =======
[mouse_expression, ~, mouse_gene_info, ~, ~, ~] = load_expression_and_regions(dataset_name,[]);

[mouse_expression, mouse_gross_region_vec, mouse_gene_info, mouse_samples2subjects, mouse_gross_region_data] = load_expression_and_regions('mouse','standard+cortex');

homologous_genes_human6_zapala(human_gene_info, mouse_gene_info)