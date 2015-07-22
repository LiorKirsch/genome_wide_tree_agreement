function [scores, symbols_out, entrez_out,label] = get_subset_scores(...
    subset_name, human6Results, human_gene_info, human6RandomResults)

splits = textscan(subset_name,'%s %d','Delimiter','-');
subset_name = splits{1}{1};
if ~isempty(splits{2})
    extra_parm = splits{2}(1);
end

switch subset_name
    case 'serotonin'
        gene_sym = {'HTR1A','HTR1B','HTR1D','HTR1E','HTR1F','HTR2A','HTR2B','HTR2C','HTR3A','HTR3B','HTR3C','HTR3D','HTR3E','HTR4','HTR5A','HTR6','HTR7','TPH1','TPH2', 'SLC6A4'};
        gene_entrez = nan(size(gene_sym));
        label = 'Serotonin';
    case 'dopamin' 
        gene_sym = {'COMT','DRD1','DRD2','DRD3','DRD4','DRD5','SLC6A3','TH'};
        gene_entrez = nan(size(gene_sym));
        label = 'Dopamin';
    case 'dopaminAndSertonin' 
        gene_sym = {'DDC', 'MAOA','MAOB','SLC18A1','SLC18A2'};
        gene_entrez = nan(size(gene_sym));
        label = 'Dopamin and Sertonin';
    case 'SOX2'
        gene_sym = {'SOX2';};
        gene_entrez = nan(size(gene_sym));
        label = 'Sox2';
    case 'PAX'
        gene_sym = {'Pax1','Pax2','Pax3','Pax4','Pax5','Pax6','Pax7','Pax8';};
        gene_entrez = nan(size(gene_sym));
        label = 'Pax';
    case 'HOX'
        gene_sym = {'Hoxa1';'Hoxa10';'Hoxa11';'Hoxa2';'Hoxa3';'Hoxa4';'Hoxa5';'Hoxa6';'Hoxa7';'Hoxa9';'Hoxb1';'Hoxb13';'Hoxb3';'Hoxb4';'Hoxb5';'Hoxb6';'Hoxb9';'Hoxc10';'Hoxc12';'Hoxc13';'Hoxc4';'Hoxc5';'Hoxc8';'Hoxc9';'Hoxd1';'Hoxd12';'Hoxd13';'Hoxd3';'Hoxd4';'Hoxd8';'Hoxd8';'Hoxd9';};
        gene_entrez = nan(size(gene_sym));
        label = 'Hox';
    case 'axon_guidance'
        gene_sym = textread('/cortex/data/gene_sets/human_axon_guidance','%s');
        gene_entrez = nan(size(gene_sym));
        label = 'Axon guidance';
    case 'cahoy_astro'
        cahoy_human = load('cahoy_human_homolgues.mat','astro_human_entrez','astro_human_hsbc');
        gene_entrez = cahoy_human.astro_human_entrez;
        gene_sym = cahoy_human.astro_human_hsbc;
        label = 'Astrocytes';
    case 'cahoy_oligo'
        cahoy_human = load('cahoy_human_homolgues.mat','oligo_human_entrez','oligo_human_hsbc');
        gene_entrez = cahoy_human.oligo_human_entrez;
        gene_sym = cahoy_human.oligo_human_hsbc;
        label = 'Oligodendrocytes';
    case 'cahoy_neuro'
        cahoy_human = load('cahoy_human_homolgues.mat','neurons_human_entrez','neurons_human_hsbc');
        gene_entrez = cahoy_human.neurons_human_entrez;
        gene_sym = cahoy_human.neurons_human_hsbc;
        label = 'Neurons';
    case 'housekeeping'
        [housekeeping.ensembl_gene_id,housekeeping.transcript_gene_id ,housekeeping.entrez ,housekeeping.symbol ] = textread('/cortex/data/gene_sets/house_keeping_genes/Erez_2003.txt','%s %s %d %s','headerlines',1);
        gene_entrez = housekeeping.entrez;
        gene_sym = housekeeping.symbol;
        label = 'Housekeeping';
    case 'age'
        load('gene_ages.mat');
        num_age = extra_parm;
        gene_sym = selectedProbesData.gene_symbols(logical(ages(:,num_age)));
        gene_entrez = selectedProbesData.entrez_ids(logical(ages(:,num_age)));
        label = agesDescription{num_age};
     
    otherwise
        error('unkown subset %s', subset_name)
end

gene_sym = upper(gene_sym);

% [~, indInList] = ismember(gene_sym,human_gene_info.gene_symbols);
% intersect_symbol = human_gene_info.gene_symbols(indInList);
% intersect_entrez = human_gene_info.entrez_ids(indInList);
% [scores, symbols_out, entrez_out] = addScoresToSubset(human6Results, human_gene_info, intersect_entrez, intersect_symbol);
[scores, symbols_out, entrez_out] = addScoresToSubset(human6Results, human_gene_info, gene_entrez, upper(gene_sym));

empirical_pvalue = getEmpiricalPvalues(scores, human6RandomResults);
empirical_pvalue = mafdr(empirical_pvalue, 'BHFDR', true);
fprintf('%s: %d / %d  (%d%%)\n', label, sum(empirical_pvalue < 0.01) , length(empirical_pvalue) ,floor(sum(empirical_pvalue < 0.01) / length(empirical_pvalue)*100) );
compareDistributions(scores, human6Results,  label, 'all human6');
end