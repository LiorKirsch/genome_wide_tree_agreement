function anaylizeTopGenes()
% analyzes the top BRO genes. (Enrichment, GO-cat...)
% 

init;
% data_results = load('results/human6GrossRegions-1-30000.mat');
data_results = load('results/human6AllRegions-1-30000.mat');
human6Results = data_results.results;
human6RandomResults = data_results.randomResults;
human_gene_info = data_results.gene_info;

% data_results = load('results/zapalaMouse-1-30000.mat');
data_results = load('results/kangAllRegions-1-30000.mat');
second_dataset_results = data_results.results;
second_dataset_random_results = data_results.randomResults;
second_dataset_genes_info = data_results.gene_info;


[ind_for_A, ind_for_B] = mapGenes(human_gene_info.gene_symbols, second_dataset_genes_info.gene_symbols, human_gene_info.entrez_ids, second_dataset_genes_info.entrez_ids);

scores_A = human6Results(ind_for_A);
scores_B = second_dataset_results(ind_for_B);
sym = human_gene_info.gene_symbols(ind_for_A);
entrez = human_gene_info.entrez_ids(ind_for_A);
% sym_B = second_dataset_genes_info.gene_symbols(ind_for_B);

mult_score = scores_A .* scores_B;
[sorted_score,sort_ind] = sort(mult_score, 'descend');
sorted_symb = sym(sort_ind);
sorted_entrez = entrez(sort_ind);
sorted_scores_A = scores_A(sort_ind);
sorted_scores_B = scores_B(sort_ind);

sorted_entrez = arrayfun(@(x) sprintf('%d',x), ...
    sorted_entrez,'UniformOutput',false);

top = 20;
for i = 1:top
    url = sprintf('http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', sorted_symb{i} );
    fprintf('%4.4g \t%s \t %s\n', sorted_score(i) , sorted_symb{i}, url );
end

cell2_ind = print_BRO_for_go_cat('cell_2_cell', sorted_symb, sorted_score);
synap_ind = print_BRO_for_go_cat('synaptic_trans', sorted_symb, sorted_score);
neu_diff_ind = print_BRO_for_go_cat('neuron differentiation', sorted_symb, sorted_score);
neu_proj_ind = print_BRO_for_go_cat('neuron projection development', sorted_symb, sorted_score);
gener_neu_ind = print_BRO_for_go_cat('generation of neurons', sorted_symb, sorted_score);
axon_ind = print_BRO_for_go_cat('axon development', sorted_symb, sorted_score);
neuro_dev_ind = print_BRO_for_go_cat('neuron development', sorted_symb, sorted_score);
substance_related = print_BRO_for_go_cat('substance related disorders', sorted_symb, sorted_score);
autistic = print_BRO_for_go_cat('autistic disorder', sorted_symb, sorted_score);
seizures = print_BRO_for_go_cat('seizures', sorted_symb, sorted_score);
epilepsy = print_BRO_for_go_cat('epilepsy', sorted_symb, sorted_score);
schizophrenia = print_BRO_for_go_cat('schizophrenia', sorted_symb, sorted_score);

create_csv({sorted_symb,sorted_entrez}, {sorted_score,sorted_scores_A,sorted_scores_B},...
    {cell2_ind,synap_ind,neu_diff_ind,neu_proj_ind,gener_neu_ind...
       axon_ind, neuro_dev_ind,schizophrenia,autistic,...
       seizures,epilepsy,substance_related}, ...
    {'symbol','entrez','BRO-score (product)','BRO-score (ABA6-2013)','BRO-score (Kang-2011)',...
    'cell-cell signaling','synaptic transmission','neuron differentiation'...
    'neuron projection development','generation of neurons','axon development','neuron development'...
    'schizophrenia','autistic disorder','seizures','epilepsy','substance related disorders'});


toppgene_functional_enrichment(sorted_symb(1:1000),'HGNC'); 
toppgene_functional_enrichment(sorted_entrez(1:1000),'ENTREZ');


drawGeneAcrossRegions(sorted_symb(1:10));

[~,sort_ind_6] = sort(human6Results,'descend');
sort_6_symb = human_gene_info.gene_symbols(sort_ind_6);
drawGeneAcrossRegions(sort_6_symb(1:10));

end

function create_csv(symb, scores, masks, headers)


    fid = fopen('bro-scores-out.csv','wb');
    print_line = '';
    for i = 1:length(headers)
        print_line = sprintf('%s%s,', print_line, headers{i} );
    end
    fprintf(fid,'%s\n', print_line(1:end-1));
    
    num_genes = length(symb{1});
    for i = 1:num_genes
        print_line = '';
        for j = 1:length(symb)
           print_line = sprintf('%s%s,', print_line, symb{j}{i} );
        end
        for j = 1:length(scores)
           print_line = sprintf('%s%g,', print_line, scores{j}(i) );
        end
        for j = 1:length(masks)
           if masks{j}(i)
                print_line = sprintf('%sX,', print_line );
           else
                print_line = sprintf('%s,', print_line);
           end
        end
        
        fprintf(fid,'%s\n', print_line(1:end-1));
    end
    fclose(fid);
end