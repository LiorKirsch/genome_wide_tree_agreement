function [scoresA, gene_symbols, entrez_ids] = addScoresToSubset(scoresA, geneInfo, subset_entrez, subset_names)
    [ind_for_A, ind_for_B] = mapGenes(geneInfo.gene_symbols,subset_names, geneInfo.entrez_ids, subset_entrez);
        
    scoresA = scoresA(ind_for_A);
    gene_symbols = geneInfo.gene_symbols(ind_for_A);
    entrez_ids = geneInfo.entrez_ids(ind_for_A);
end