function [ind_for_A, ind_for_B] = mapGenes(gene_symbols_A, gene_symbols_B, entrez_A, entrez_B)
% return the indices of the genes that appear in both A and B
% checks for the same gene name and adds genes with the same entrez

   gene_symbols_A = upper(gene_symbols_A);
   gene_symbols_B = upper(gene_symbols_B);

   [intersect_symbols, gene_symbols_ind_for_A, gene_symbols_ind_for_B] = intersect(gene_symbols_A, gene_symbols_B);
   ismember_A_genes = ismember(gene_symbols_A, intersect_symbols);
   ismember_B_genes = ismember(gene_symbols_B, intersect_symbols);
   intersect_entrez = intersect(entrez_A, entrez_B);
   ismember_B_entrez = ismember(entrez_B, intersect_entrez);
   ismember_A_entrez = ismember(entrez_A, intersect_entrez);

    
   gene_in_entrez_not_in_symbols_ind_A = ~ismember_A_genes & ismember_A_entrez;
   gene_in_entrez_not_in_symbols_ind_B = ~ismember_B_genes & ismember_B_entrez;
   
   extra_entrez = intersect( entrez_A(gene_in_entrez_not_in_symbols_ind_A) , entrez_B(gene_in_entrez_not_in_symbols_ind_B) );
   
%    gene_in_entrez_not_in_symbols = entrez_A(gene_in_entrez_not_in_symbols_ind);
   [~, extra_ind_entrez_A] = ismember(extra_entrez, entrez_A);
   [~, extra_ind_entrez_B] = ismember(extra_entrez, entrez_B);
   
   ind_for_A = [gene_symbols_ind_for_A; extra_ind_entrez_A];
   ind_for_B = [gene_symbols_ind_for_B; extra_ind_entrez_B];
   
end