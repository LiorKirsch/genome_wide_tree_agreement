function checkDevlopmentAgreement()


 kangAges = {'4-8pcw', '8-10pcw', '10-13pcw', '13-16pcw', '16-19pcw', '19-24pcw', '24-38pcw', '0-6mo', '6-12mo', '1-6y', '6-12y', '12-20y', '20-40y', '40-60y', '60+y'};

 fileName = 'results/simpleMeasurementKangAges.mat';
 kangAgesAgreement = load(fileName);
 
 gene_name = 'HOXA1';
 indices = strcmp(kangAgesAgreement.developing_genes_info.gene_symbols , gene_name);
 gene_bro_score = nan(length(kangAges),1);
    for i =4:length(kangAges)
        results = kangAgesAgreement.brainspanResults{i};
        gene_bro_score(i) = results(indices);

        
    end
 
figure('name',sprintf('Bro scores during dev - %s', gene_name) );
plot(gene_bro_score);
xticklabel_rotate(1:length(kangAges),45,kangAges);

   

hoxGenesSymbols = {'Hoxa1';'Hoxa10';'Hoxa11';'Hoxa2';'Hoxa3';'Hoxa4';'Hoxa5';'Hoxa6';'Hoxa7';'Hoxa9';'Hoxb1';'Hoxb13';'Hoxb3';'Hoxb4';'Hoxb5';'Hoxb6';'Hoxb9';'Hoxc10';'Hoxc12';'Hoxc13';'Hoxc4';'Hoxc5';'Hoxc8';'Hoxc9';'Hoxd1';'Hoxd12';'Hoxd13';'Hoxd3';'Hoxd4';'Hoxd8';'Hoxd8';'Hoxd9';};

figure('name', 'Hox family Bro scores during dev' );
showGroupOfGenesByAge(kangAges, kangAgesAgreement,hoxGenesSymbols);

end

function showGroupOfGenesByAge(kangAges, kangAgesAgreement,intrestGroupGeneNames)

    [~, indInList] = ismember(upper(intrestGroupGeneNames),kangAgesAgreement.developing_genes_info.gene_symbols );

    indInList = indInList( indInList > 0);
    indInList  = unique(indInList);
    
    gene_bro_score = nan(length(kangAges),length(indInList) );
    for i =4:length(kangAges)
        results = kangAgesAgreement.brainspanResults{i};
        gene_bro_score(i,:) = results(indInList);


    end

    plot(1:length(kangAges), gene_bro_score,'.');

end