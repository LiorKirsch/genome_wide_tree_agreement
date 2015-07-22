% Analyzes the relations between the BRO-score and gene-age
%

init;
set(gca,'FontSize', 22);

load('data_matfile/all_ages.mat');

% agesDes ={'Life before LCA of Cellular_organisms - Cellular organisms';'Cellular organisms - Eukaryota';'Eukaryota - Opisthokonta';'Opisthokonta - Holozoa';'Holozoa - Metazoa';'Metazoa - Eumetazoa';'Eumetazoa - Bilateria';'Bilateria - Deuterostomia';'Deuterostomia - Chordata';'Chordata - Olfactores';'Olfactores - Craniata';'Craniata - Euteleostomi';'Euteleostomi - Tetrapoda';'Tetrapoda - Amniota';'Amniota - Mammalia';'Mammalia - Eutheria';'Eutheria - Boreoeutheria';'Boreoeutheria - Euarchontoglires';'Euarchontoglires - Primates'};
num_age = find(strcmp('Opisthokonta - Holozoa', agesDescription));
gene_sym = selectedProbesData.gene_symbols(logical(ages(:,num_age)));
gene_entrez = selectedProbesData.entrez_ids(logical(ages(:,num_age)));

[scores, symbols_out, entrez_out] = addScoresToSubset(human6Results, human_gene_info, gene_entrez, upper(gene_sym));



[scores , sort_ind] = sort(scores,'descend');

disp( symbols_out(sort_ind(1:5)));
disp( scores(1:5));


gene_sym = selectedProbesData.gene_symbols(sum(ages,2)==1);
gene_entrez = selectedProbesData.entrez_ids(sum(ages,2)==1);
age_for_genes = ages(sum(ages,2)==1,:)*((1:19)');

[ind_for_A, ind_for_B] = mapGenes(human_gene_info.gene_symbols,upper(gene_sym), human_gene_info.entrez_ids, gene_entrez);
score_per_gene = human6Results(ind_for_A);
age_per_gene = age_for_genes(ind_for_B);

figure;
% hold on;
% mean_per_age = accumarray(age_per_gene, score_per_gene, [], @mean);
% std_per_age = accumarray(age_per_gene, score_per_gene, [], @std);
% % H=shadedErrorBar(1:length(mean_per_age),mean_per_age,std_per_age,'-b') ; %,lineProps,transparent)
% plot(mean_per_age)
% % errorbar(1:length(mean_per_age),mean_per_age,std_per_age);

median_per_age = accumarray(age_per_gene, score_per_gene, [], @median);
ploth = plot(median_per_age,'k');
set(ploth,'LineWidth',  2 );
set(ploth,'Color',  [1.0000,0.2,0] );

% legend('mean BRO score','median BRO score' , 'location','best');
% legend('boxoff');
ylabel('Median BRO agreement'); xlabel('Gene age');

x_marks = cell(length(agesDescription) +1,1);
for i = 1:length(agesDescription)
    split = strsplit(agesDescription{i} , '-');
    x_marks{i} = strtrim(split{1});
    x_marks{i+1} = strtrim(split{2});
end

x_marks{1} =  '';
    ax = gca;
    ax.XTick = (1:length(x_marks)) - 0.5;
    ax.XTickLabel = x_marks;
    ax.XTickLabelRotation	=45;
    
set(gca,'Fontsize',16);         %set(gca,'xscale','log');    
%% 
    
data_results = load('results/human6AllRegions-1-30000.mat');
human6Results = data_results.results;
human6RandomResults = data_results.randomResults;
human_gene_info = data_results.gene_info;
all_samples_color = [0.8, 0.8,0.8];
xlimit = [-0.13,0.8]; % full tree
ylimit = [0, 0.25] ; % full tree


age_scores = {}; age_labels = {};
for i =[1,19]
    [age_scores_i,~,~,age_labels_i] = get_subset_scores(...
    sprintf('age-%d',i ), human6Results, human_gene_info,human6RandomResults);
    age_scores = cat(2,age_scores, {age_scores_i});
    age_labels = cat(2,age_labels, {age_labels_i});
end
age_colors = autumn(length(age_labels) +1);
age_colors = age_colors(1:end-1,:);

 
age_labels = cellfun(@(x) strsplit(x,' - '), age_labels, 'uniformOutput',false);
age_labels = cellfun(@(x) x{2}, age_labels, 'uniformOutput',false);


draw_groups_dist('Age',[{human6Results}, age_scores], ...
    [{'All'}, age_labels], [all_samples_color;age_colors], xlimit,ylimit,true);

compareDistributions(age_scores{1}, age_scores{2},  age_labels{1}, age_labels{2});

