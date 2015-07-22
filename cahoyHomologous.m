% Creates The human marker list for the Cahoy neuronal markers dataset
%

addpath('/home/lab/lior/Projects/general use functions');

% load the human entrez from the allen 6 human dataset
load('/cortex/data/microarray/human/Hawrylycz_2012/easyFormatHumanData.mat','selectedProbesData');
humanEntrez = selectedProbesData.entrez_ids;

% load the list of genes from cahoy above a threshold of x Fold
fold_threshold = 10;
[neurons_names, neurons_entrez, astro_names, astro_entrez, oligo_names, oligo_entrez] =  load_cahoy_above_fold(fold_threshold);


% get the human homologues of the cahoy mouse genes
[~, ~, ~, neurons_human_hsbc, neurons_human_entrez] = findHomolgousMouseHuman(neurons_entrez, humanEntrez);
[~, ~, ~, astro_human_hsbc, astro_human_entrez] = findHomolgousMouseHuman(astro_entrez, humanEntrez);
[~, ~, ~, oligo_human_hsbc, oligo_human_entrez] = findHomolgousMouseHuman(oligo_entrez, humanEntrez);

save('data_matfile/cahoy_human_homologous.mat', 'neurons_human_hsbc', 'neurons_human_entrez', 'astro_human_hsbc', 'astro_human_entrez', 'oligo_human_hsbc', 'oligo_human_entrez','fold_threshold');