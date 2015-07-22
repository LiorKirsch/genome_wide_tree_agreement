function cat_mask = print_BRO_for_go_cat(go_cat, gene_symbols, gene_scores)

switch go_cat
    case 'schizophrenia'
        filename = 'go_cat/schizophrenia.csv';
    case 'epilepsy'
        filename = 'go_cat/epilepsy.csv';
    case 'autistic disorder'
        filename = 'go_cat/autistic disorder.csv';
    case 'seizures'
        filename = 'go_cat/seizures.csv';
    case 'substance related disorders'
        filename = 'go_cat/substance related disorders.csv';
    case 'synaptic_trans'
        filename = 'go_cat/synaptic_trans.csv';
    case 'cell_2_cell'
        filename = 'go_cat/cell_2_cell.csv';
    case 'neuron differentiation'
        filename = 'go_cat/neuron differentiation.csv';
    case 'neuron projection development'
        filename = 'go_cat/neuron projection development.csv';
    case 'generation of neurons'
        filename = 'go_cat/generation of neurons.csv';
    case 'axon development'
        filename = 'go_cat/axon development.csv';
    case 'neuron development'
        filename = 'go_cat/neuron development.csv';
    otherwise
        error('unkown category %s', go_cat);
end

z = textscan(fopen(filename),'%q %q %q','Delimiter',',','HeaderLines',1);
cat_entrez = z{1};
cat_symb = z{2};
cat_description = z{3};


cat_mask = ismember(gene_symbols, cat_symb);
[~,inds_a] = intersect(gene_symbols, cat_symb);

intersect_symb = gene_symbols(inds_a);
intersect_scores = gene_scores(inds_a);
[intersect_scores,sort_inds] = sort(intersect_scores,'descend');
intersect_symb = intersect_symb(sort_inds);

fid = fopen(sprintf('go_cat/%s_out.txt', go_cat) ,'wb');
for i = 1:length(intersect_scores)
     url = sprintf('http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', intersect_symb{i} );

    fprintf(fid, '%g\t%s\t%s\n', intersect_scores(i), intersect_symb{i},url);
end

