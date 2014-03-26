function mapMouseToHuman()
addpath('/home/lab/gal/develop/matlab');

    x = load('/home/lab/gal/Projects/Limor/Data/ABA_expression_data.mat');
    [mapping, genes_not_found] = createAllenToEntrezMapping(x.name_of_genes,'allen');
    
    emptyCells = cellfun(@isempty,mapping(:,1) );
    mapping(emptyCells,1) = {nan};
    
    allen_mouse_entrez = mapping(:,1);
    allen_mouse_entrez = cell2mat(allen_mouse_entrez);
%     genesWithEntrez = ~isnan(allen_mouse_entrez);
%     allen_mouse_entrez = allen_mouse_entrez( genesWithEntrez );
    fid = fopen('allen_mouse_entrez.txt','w');
    fprintf(fid,'%d\n',allen_mouse_entrez);
    fclose(fid);
    

    % -------------------------------------------
    % run R code to map the mouse entrez to human 
    % 
    % do find and replace  Na with  NaN
    %
    % -------------------------------------------
    
    [human_housekeeping.Ensembl_Gene_ID,	human_housekeeping.Ensembl_Transcript_ID,human_housekeeping.EntrezGene_ID, human_housekeeping.HGNC_symbol] = textread('/cortex/data/gene_sets/house_keeping_genes/Erez_2003.txt','%s %s %d %s','delimiter','\t','headerlines',1);
    [~, mouse_mgi, mouse_wikiname, mouse_entrez, human_hsbc, human_entrez] = textread('Mouse2Human.data', '%s %q %q %d %q %d','headerlines',1);
    
   
        
    h = load('/home/lab/lior/Projects/individual variability/easyFormatHumanData.mat','selectedProbesData');
    
   
    [ind_for_human_data, ind_for_mouse_data, ~,relevent_genes_with_orthologs] = getTheJointIndcies(human_entrez,mouse_entrez, h.selectedProbesData.entrez_ids, allen_mouse_entrez);
    numberOfJointGenes = sum(relevent_genes_with_orthologs);
    
%     test1 = h.selectedProbesData.gene_symbols(ind_for_human_data);
%     test2 = human_hsbc(relevent_genes_with_orthologs);
%     testB1 = x.name_of_genes(ind_for_mouse_data)';
%     testB2 = mouse_wikiname(relevent_genes_with_orthologs);

  
% get the expression and the region vec from the dataset then slice and sort it by the
% order that is common to both mouse and human.

    [human_expression, human_gross_region_vec, human_gene_info, human_samples2subjects, human_gross_region_names] = load_expression_and_regions('human6',[]);
    human_genes_symbols = human_gene_info.gene_symbols;
    human_expression = human_expression(:,ind_for_human_data);
    human_genes_symbols = human_genes_symbols(ind_for_human_data);

    [mouse_expression, mouse_gross_region_vec, mouse_genes_info, ~, mouse_gross_region_names] = load_expression_and_regions('mouse',[]);
    mouse_genes_symbols = mouse_genes_info.gene_symbols;
    mouse_expression = mouse_expression(:,ind_for_mouse_data);
    mouse_genes_symbols = mouse_genes_symbols(ind_for_mouse_data)';

%     genesArray = {'Adcy2','Neud4','Gda','Chn1','Sprn','Kctd16','Kcnf1','Mmp17','Ak5','Fhl2','Lmo4','Csmd1','Gria2','Bai2' };
%     genesArray = {'Zfp521','Fnbp1','E130112L23Rik','AF529169','Serpinb1a','Elovl5','Osbpl9','Akap12','Syt2','Clcn5','Rec8L1','Acyp2','Agt','Dnm3','Mcam','Ece2','Mtss1','S100b','Phka1','Zfp423'};
%     genesArray = {'Hod','Serpinf1','Slc12a8','Steap2','Chrm2','Col5a3','Tm2d3','Echs1','2410022L05Rik','Kirrel2','Kdelr3','Kif5a','Cdkl1'};
% 
%     for i = 1:length(genesArray)
%         geneName =  genesArray{i};
%         subplot(1,2,1);
%         drawGeneExpressionInGrossRegions(geneName, mouse_genes_symbols, human_expression, human_gross_region_vec, human_gross_region_names);
%         title([geneName , ' human expression'],'fontsize',20); subplot(1,2,2);
%         drawGeneExpressionInGrossRegions(geneName, mouse_genes_symbols, mouse_expression, mouse_gross_region_vec, mouse_gross_region_names);
%         title([geneName , ' mouse expression'],'fontsize',20);
%         saveFigure(gcf, [geneName, '.png'], 'png')
%     end
    
    human_scores_all = nan(6,size(human_expression,2));
    for i=1:size(human_samples2subjects,2)
        current_human_samples = human_samples2subjects(:,i);
        human_scores_all(i,:) = scatterMouseHumanCorr(human_expression(current_human_samples,:), human_gross_region_vec(current_human_samples), mouse_expression,mouse_gross_region_vec, human_genes_symbols,mouse_genes_symbols);
        title(sprintf('human %d (%d samples)' , i, sum(current_human_samples) ),'fontsize',20);
    end
    
    all_human_scores = corr(human_expression,human_gross_region_vec,'type','Spearman');
    meanHumanScore = mean(human_scores_all,1)';
    stdHumanScore = std(human_scores_all,0,1)';
    mouse_scores = corr(mouse_expression, mouse_gross_region_vec,'type','Spearman');
    notnan = ~isnan(mouse_scores);
    mouse_scores = mouse_scores(notnan);
    all_human_scores = all_human_scores(notnan);
    meanHumanScore = meanHumanScore(notnan);
    human_genes_symbols_notnan = human_genes_symbols(notnan);
    mouse_genes_symbols_notnan = mouse_genes_symbols(notnan);
    
    generateCombCorrelationScore(all_human_scores, mouse_scores, human_genes_symbols_notnan, 'mouse_humans');
    
    stdStrings = arrayfun(@num2str, stdHumanScore , 'unif', 0);
    stdStrings = stdStrings(notnan);
    symbolsWithStd = strcat( human_genes_symbols_notnan,  repmat({'  \pm '}, size(human_genes_symbols_notnan)) ,stdStrings);
     scatterWithString(meanHumanScore, mouse_scores, symbolsWithStd, 'human correlation scores', 'mouse correlation scores', 'spearman correlation with the neural tube ordering');
     drawContours();
     
    results_low_low = getPercintileGenes(meanHumanScore, mouse_scores, [0 0.2], [0 0.2]);
    results_high_high = getPercintileGenes(meanHumanScore, mouse_scores, [0.8 1], [0.8 1]);
    results_high_low = getPercintileGenes(meanHumanScore, mouse_scores, [0.8 1], [0 0.2]);
    
    writeCellToFile('human__low_human_low_mouse.txt', human_genes_symbols_notnan(results_low_low) );
    writeCellToFile('human__high_human_high_mouse.txt', human_genes_symbols_notnan(results_high_high) );
    writeCellToFile('human__high_human_low_mouse.txt', human_genes_symbols_notnan(results_high_low) );
    writeCellToFile('human__background.txt', human_genes_symbols_notnan );
    
    
    writeCellToFile('mouse__low_human_low_mouse.txt', mouse_genes_symbols_notnan(results_low_low) );
    writeCellToFile('mouse__high_human_high_mouse.txt', mouse_genes_symbols_notnan(results_high_high) );
    writeCellToFile('mouse__high_human_low_mouse.txt', mouse_genes_symbols_notnan(results_high_low) );
    writeCellToFile('mouse__background.txt', mouse_genes_symbols_notnan );
    
    [~,sortByMult] = sort(mouse_scores.* meanHumanScore);
    writeCellToFile('mouse__sortByMult.txt', mouse_genes_symbols_notnan(sortByMult) );
    writeCellToFile('human__sortByMult.txt', human_genes_symbols_notnan(sortByMult) );
end

function drawGeneExpressionInGrossRegions(geneName, genes_symbols, experimentsDataMatrix, gross_region_vec, regionNames)
    geneIndex = strcmp(genes_symbols, geneName);
    geneExpression = experimentsDataMatrix(:,geneIndex);
    numberOfRegions = max(gross_region_vec);
    for i=1:numberOfRegions
       meanExpression(i) =  mean(geneExpression(gross_region_vec == i ));
       stdExpression(i) = std(geneExpression(  gross_region_vec == i  ));
    end
    regionColors = repmat([0 0 1],numberOfRegions,1);
    
%     drawBars(meanExpression,stdExpression,regionNames,regionColors,'mean expression');
    
    errorbar(1:length(meanExpression), meanExpression, stdExpression);
    xticklabel_rotate(1:length(meanExpression),40,regionNames, 'fontSize',14) ;
    ylabel('mean expression','fontSize',14);
end

function results = getPercintileGenes(human_scores, mouse_scores, human_range, mouse_range)

    [~,h_sort_ind] = sort(human_scores);
    [~,reverse_h_sort_ind] = sort(h_sort_ind);
    [~,m_sort_ind] = sort(mouse_scores);
    [~,reverse_m_sort_ind] = sort(m_sort_ind);
    results_h = human_range(1) * length(human_scores) <= reverse_h_sort_ind & reverse_h_sort_ind <= human_range(2)* length(human_scores);
    results_m = mouse_range(1) * length(mouse_scores) <= reverse_m_sort_ind & reverse_m_sort_ind <= mouse_range(2)* length(mouse_scores);
    results = results_h & results_m ;
end

function generateCombCorrelationScore(corr1, corr2, strings, file_name)
    positive = (corr1 > 0) & (corr2 > 0);
    negative = (corr1 < 0) & (corr2 < 0);
    
    scores = corr1(positive) .* corr2(positive);
    [sortedScores, ind] = sort(scores,'descend');
    sortedStrings = strings(positive);
    sortedStrings = sortedStrings(ind);
    fid = fopen( sprintf('positive_%s.txt', file_name),'w');
    for i = 1:length(sortedScores)
       fprintf( fid, '%g, %s\n', sortedScores(i), sortedStrings{i} );
    end
    fclose(fid);
    
    scores = corr1(negative) .* corr2(negative);
    [sortedScores, ind] = sort(scores,'descend');
    sortedStrings = strings(negative);
    sortedStrings = sortedStrings(ind);
    fid = fopen( sprintf('negative_%s.txt', file_name),'w');
    for i = 1:length(sortedScores)
       fprintf( fid, '%g, %s\n', sortedScores(i), sortedStrings{i} );
    end
    fclose(fid);
    
end

function human_scores = scatterMouseHumanCorr(human_expression, human_gross_region_vec, mouse_expression,mouse_gross_region_vec,human_genes_symbols,mouse_genes_symbols)
        human_scores = corr(human_expression, human_gross_region_vec,'type','Spearman');
        
        human_random_scores = getRandomCorrelations(human_expression,human_gross_region_vec,5);
    
        mouse_scores = corr(mouse_expression, mouse_gross_region_vec,'type','Spearman');

        isNaN = isnan(human_scores) | isnan(mouse_scores);


        mouse_random_scores = getRandomCorrelations(mouse_expression,mouse_gross_region_vec,5);

        isNaN = isnan(human_random_scores) | isnan(mouse_random_scores);
        randCov = cov([human_random_scores(~isNaN), mouse_random_scores(~isNaN)]);
        randMean = mean([human_random_scores(~isNaN), mouse_random_scores(~isNaN)],1);

        figure;
        hold on;
        scatter(human_scores, mouse_scores,10, [ .5 .5 .5], 'filled');
        error_ellipse(randCov,randMean);
        error_ellipse(2*randCov,randMean,'style','--');
        xlabel('human spearman','fontsize',20);
        ylabel('mouse spearman','fontsize',20);

        hoxGenes = {'Hoxa1';'Hoxa10';'Hoxa11';'Hoxa2';'Hoxa3';'Hoxa4';'Hoxa5';'Hoxa6';'Hoxa7';'Hoxa9';'Hoxb1';'Hoxb13';'Hoxb3';'Hoxb4';'Hoxb5';'Hoxb6';'Hoxb9';'Hoxc10';'Hoxc12';'Hoxc13';'Hoxc4';'Hoxc5';'Hoxc8';'Hoxc9';'Hoxd1';'Hoxd12';'Hoxd13';'Hoxd3';'Hoxd4';'Hoxd8';'Hoxd8';'Hoxd9';};
        scatterSelectedList(human_scores, mouse_scores,mouse_genes_symbols, hoxGenes, [ 1 0 0]);

        % could not find gene Otx2
        scatterSelectedList(human_scores, mouse_scores,mouse_genes_symbols, {'Gbx2'}, [ 0 1 0]);

        human_axon_gudiance_genes = {'CDK5  ';'HRAS  ';'GSK3B  ';'CXCR4  ';'PPP3C';'RASA1';'ERK1_2  ';'RAC1  ';'CDC42  ';'PAK1  ';'PAK2  ';'RHOA  ';'ROCK  ';'GNAI  ';'MET';'EPHA1';'EPHA2';'EPHA3';'EPHA4';'EPHA5';'EPHA6';'EPHA7';'EPHA8';'EPHB1';'EPHB2';'EPHB3';'EPHB4';'EPHB6  ';'EFNA  ';'EFNB  ';'FYN  ';'ITGB1  ';'PTK2';'PAK3';'PAK4  ';'PAK6  ';'PAK7  ';'LIMK1  ';'LIMK2  ';'CFL  ';'PPP3R';'SEMA4  ';'SEMA7';'L1CAM  ';'PLXNC  ';'ABL1  ';'NRP1  ';'ROBO1  ';'ROBO2  ';'ROBO3  ';'DCC  ';'PLXNA  ';'PLXNB  ';'SLIT1  ';'SLIT2  ';'SEMA3  ';'SEMA5  ';'SEMA6  ';'NTN1  ';'NTN3  ';'NTN4  ';'SLIT3  ';'ABLIM  ';'UNC5  ';'NTNG1  ';'NGL1  ';'RGS3  ';'NGEF';'SRGAP  ';'FPS';'DPYSL2';'DPYSL5';'RHOD  ';'RND1  ';'ARHGEF12';'KRAS';'NRAS  ';'RAC2  ';'RAC3  ';};
        scatterSelectedList(human_scores, mouse_scores,human_genes_symbols, human_axon_gudiance_genes, [ 0 0 1]);

%         [human_housekeeping.Ensembl_Gene_ID,	human_housekeeping.Ensembl_Transcript_ID,human_housekeeping.EntrezGene_ID, human_housekeeping.HGNC_symbol] = textread('/cortex/data/gene_sets/house_keeping_genes/Erez_2003.txt','%s %s %d %s','delimiter','\t','headerlines',1);
%         scatterSelectedList(human_scores, mouse_scores,human_genes_symbols, human_housekeeping.HGNC_symbol, [ 0.5 0.5 0]);
%        legend('genes with orthologs','cov','2*cov','hox gene family','Gbx2','axon gudiance genes','human housekeeping','Location','NorthWest');        
 
        legend('genes with orthologs','cov','2*cov','hox gene family','Gbx2','axon gudiance genes','Location','NorthWest');
        

         scatterWithString(human_scores, mouse_scores, mouse_genes_symbols, 'human correlation scores', 'mouse correlation scores', 'spearman correlation with the neural tube ordering')
    %     corrBetween = corr(relevent_human_scores(~isNaN), relevent_mouse_scores(~isNaN),'type','Spearman');
end

function getGenesOnBorders(human_scores_all, mouse_scores, human_genes_symbols, mouse_genes_symbols)
    mean_human_score = mean(human_scores_all,1);
    
end

function scatterSelectedList(human_scores, mouse_scores, fullList, selectedList, color)
    selected_genes_ind = ismember(fullList, selectedList);
    scatter(human_scores(selected_genes_ind), mouse_scores(selected_genes_ind),80, color, 'filled');
    xlim([-1 1]);
    ylim([-1 1]);
    text( human_scores(selected_genes_ind) ,mouse_scores(selected_genes_ind) , fullList(selected_genes_ind), 'horizontal','left', 'vertical','bottom','fontsize',14);
end

function drawContours()
    hold on;
    X = -1:0.01:1;
    Y = -1:0.01:1;
    [X,Y] = meshgrid(X,Y);
    Z = X.*Y;
    contour(X,Y,Z,50);
    hold off;
end
function scatterWithString(xAxis, yAxis, strings, xtitle, ytitle, overalltitle)
    figure;
    scatter(xAxis, yAxis,'.');
    xlabel(xtitle,'fontsize',20);
    ylabel(ytitle,'fontsize',20);
    title(overalltitle,'fontsize',20);
    text( xAxis ,yAxis , strings, 'horizontal','left', 'vertical','bottom','fontsize',10);

end

function randCorr = getRandomCorrelations(data,data2,repeat)

    numberOfSamples = size(data,1);
    numberOfFeatures = size(data,2);
    
    randCorr = nan(numberOfFeatures*repeat,1);
    for i = 1:repeat
        randomPermutation = randperm(size(data,1));
        randCorr( (i-1)*numberOfFeatures +1: i*numberOfFeatures )   = corr(data, data2(randomPermutation,:),'type','Spearman');
    end
    
end

function [ind_for_human_data, ind_for_mouse_data, ind_for_mapping,relevent_genes_with_orthologs] = getTheJointIndcies(human_joint_ids,mouse_joint_ids, human_ids, mouse_ids)
 % since we have many lists we will order all of them to the order of
    % the mapping
    
%     relevent_genes_with_orthologs = ismember(human_entrez, h.selectedProbesData.entrez_ids);
%     relevent_genes_with_orthologs = relevent_genes_with_orthologs & ismember(mouse_entrez, allen_mouse_entrez);
    
    jointVectorSize = length(mouse_joint_ids);
    
    [~, ind_human,~] = intersect(human_joint_ids , human_ids);
    relevent_genes_human = false(jointVectorSize,1); relevent_genes_human(ind_human) =true;
    
    [~, ind_mouse,~] = intersect(mouse_joint_ids, mouse_ids);
    relevent_genes_mouse = false(jointVectorSize,1); relevent_genes_mouse(ind_mouse) =true;
    
    relevent_genes_with_orthologs = relevent_genes_human & relevent_genes_mouse;
    
    [human_relevant_entrez, ind_for_human_data, ind_for_mapping]= intersect_sort_by_B(human_ids, human_joint_ids(relevent_genes_with_orthologs) );
    
    [mouse_relevant_entrez, ind_for_mouse_data, ind_for_mapping]= intersect_sort_by_B(mouse_ids, mouse_joint_ids(relevent_genes_with_orthologs) );
end

function cellNumberToArray(cellArray)

    outputNum  = nan(length(cellArray),1);
    
    for i = 1:length(cellArray)
       val =  cellArray{i};
       if ~strcmp(val, 'NA')
           outputNum(i)  = double(val);
       end
    end

end