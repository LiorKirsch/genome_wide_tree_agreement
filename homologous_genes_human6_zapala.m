function homologous_genes_human6_zapala(human6_gene_info, zapala_gene_info)

addpath('/cortex/code/matlab/homologous_gene_mapping/');

% zapala_entrez_id = num2cell(zapala_gene_info.entrez_ids);
% zapala_entrez_id = cellfun(@num2str, zapala_entrez_id, 'UniformOutput', false);
zapala_entrez_id = cellfun(@num2str, zapala_gene_info.entrez_ids, 'UniformOutput', false);
human6_entrez_id = num2cell(human6_gene_info.entrez_ids);
human6_entrez_id = cellfun(@num2str, human6_entrez_id, 'UniformOutput', false);
[gene_to_group_matrix_human, gene_to_group_matrix_mouse] =  gene_to_homolog_group('human','mouse_laboratory', human6_entrez_id, 'entrez_gene_ID',zapala_entrez_id,'entrez_gene_ID');


one_to_one_match = (sum(gene_to_group_matrix_mouse,1) == 1) & (sum(gene_to_group_matrix_human,1) == 1 );
m_1to1 = (1:size(gene_to_group_matrix_mouse,1)) * gene_to_group_matrix_mouse(:,one_to_one_match); 
h_1to1 = (1:size(gene_to_group_matrix_human,1)) * gene_to_group_matrix_human(:,one_to_one_match); 
m_1to1 = zapala_gene_info.gene_symbols(m_1to1);
h_1to1 = human6_gene_info.gene_symbols(h_1to1);

% find cells with only 1 element
one2one_match = cellfun(@(x) numel(x) ==1 , homologous_genes_symbols);

entrez_with_homolog = homologous_genes_entrez(one2one_match);
entrez_with_homolog = cellfun(@(x) x{1} , entrez_with_homolog, 'UniformOutput', false);

symbols_with_homolog = homologous_genes_symbols(one2one_match);
symbols_with_homolog = cellfun(@(x) x{1} , symbols_with_homolog, 'UniformOutput', false);

c = ismember(zapala_entrez_id, entrez_with_homolog);
d = ismember(zapala_gene_info.gene_symbols, symbols_with_homolog);
b1 = human6_gene_info.gene_symbols(one2one_match);
b2 = zapala_gene_info.gene_symbols(c);

end