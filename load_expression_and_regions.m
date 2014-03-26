function [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names, physicalLocation] = load_expression_and_regions(dataname, parms)

physicalLocation = [];

 switch dataname
      case 'human6', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names,physicalLocation] = load_human6(parms);
      case 'human6AllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names,physicalLocation] = load_human6('all regions');
     case 'human6Cortex', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names,physicalLocation] = load_human6_cortex(parms);
     case 'human6CortexAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names,physicalLocation] = load_human6_cortex('all regions');
     case 'kang', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_kang(parms);
     case 'kangAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_kang_all_regions(parms,'all');
    case 'kangCortex', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_kang_cortex(parms);
    case 'kangCortexAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_kang_all_regions(parms,'cortex');
      case 'mouse', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_mouse(parms);
      case 'brainspan',
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_brainspan(parms);
    end

end

% ===========================================================
function  [expression, gross_region_vec, gene_info,samples2subjects,gross_structures_names, physicalLocation] = load_human6( parms)
% 
    
    persistent local_experimentsDataMatrix
    persistent local_experimentsLocationMatrix    
    persistent local_humanOntology
    persistent local_gene_info
    persistent local_samples2subjects
    persistent local_grossStructures
    persistent local_physicalLocation
    

    if isempty(local_experimentsDataMatrix)
        fprintf('human6: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability','data');
        ontology_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability');
        filename = 'easyFormatHumanData.mat';
        fullname = fullfile(data_dirname, filename);
        load(fullname, 'experimentsLocationMatrix', 'experimentsDataMatrix','selectedProbesData','experimentsSubjectMatrixLogical','mri_voxel_xyz'); %#ok
        
        % Load the tree
        filename = 'humanOntologyObject.mat';
        load(fullfile(ontology_dirname, filename), 'humanOntology'); %#ok
        
        grossStructures = {'Frontal Lobe';'Cingulate gyrus';'hippocampal formation';'parahippocampal gyrus';...
                           'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';'Amygdala';'Basal Forebrain';...
                           'Striatum';'Hypothalamus';'Dorsal Thalamus';'Mesencephalon';'Cerebellar Cortex';...
                           'Pontine Tegmentum';'Myelencephalon'};

        local_humanOntology = humanOntology;
        local_experimentsDataMatrix = experimentsDataMatrix;
        local_experimentsLocationMatrix = experimentsLocationMatrix;
        local_gene_info = selectedProbesData;
        local_samples2subjects = experimentsSubjectMatrixLogical;
        local_grossStructures = grossStructures;
        local_physicalLocation = mri_voxel_xyz;
    end
    
    if strcmp(parms,'all regions')
        gross_region_vec = local_experimentsLocationMatrix * (1:size(local_experimentsLocationMatrix,2) )';
        gross_structures_names = local_humanOntology.structureLabels(:,4);
    else
        gross_region_vec = fine_to_gross(local_experimentsLocationMatrix, local_humanOntology,local_grossStructures);
        gross_structures_names = local_grossStructures;
    end
    relevant_samples = (gross_region_vec>0);        
    expression = local_experimentsDataMatrix(relevant_samples,:);
    gross_region_vec = gross_region_vec(relevant_samples);
    gene_info = local_gene_info;
    samples2subjects = local_samples2subjects(relevant_samples,:);
    physicalLocation = local_physicalLocation(relevant_samples,:);
end

% ===========================================================
function  [expression, gross_region_vec, gene_info,samples2subjects,gross_structures_names, physicalLocation] = load_human6_cortex( parms)
% 
    
    persistent local_experimentsDataMatrix
    persistent local_experimentsLocationMatrix    
    persistent local_humanOntology
    persistent local_gene_info
    persistent local_samples2subjects
    persistent local_grossStructures
    persistent local_physicalLocation
    

    if isempty(local_experimentsDataMatrix)
        fprintf('human6Cortex: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability','data');
        ontology_dirname = fullfile('/','home','lab', 'lior', 'Projects', 'individual variability');
        filename = 'easyFormatHumanData.mat';
        fullname = fullfile(data_dirname, filename);
        load(fullname, 'experimentsLocationMatrix', 'experimentsDataMatrix','selectedProbesData','experimentsSubjectMatrixLogical','mri_voxel_xyz'); %#ok
        
        % Load the tree
        filename = 'humanOntologyObject.mat';
        load(fullfile(ontology_dirname, filename), 'humanOntology'); %#ok
        
        grossStructures = {'Frontal Lobe';'Occipital Lobe';'Parietal Lobe';'Temporal Lobe';};

        local_humanOntology = humanOntology;
        local_experimentsDataMatrix = experimentsDataMatrix;
        local_experimentsLocationMatrix = experimentsLocationMatrix;
        local_gene_info = selectedProbesData;
        local_samples2subjects = experimentsSubjectMatrixLogical;
        local_grossStructures = grossStructures;
        local_physicalLocation = mri_voxel_xyz;
    end
    
    if strcmp(parms,'all regions')
        all_gross_structures_names = local_humanOntology.structureLabels(:,4);
        cortexIndices = ismember(all_gross_structures_names, local_grossStructures);
        childsParentMatrix = local_humanOntology.allChildNodes;
        all_sub_region_cortex = any( childsParentMatrix(cortexIndices,:), 1);
        experimentsLocationMatrix_cortex = local_experimentsLocationMatrix(:, all_sub_region_cortex);
        
        gross_region_vec = experimentsLocationMatrix_cortex * (1:size(experimentsLocationMatrix_cortex,2) )';
        gross_structures_names = all_gross_structures_names(all_sub_region_cortex);
    else
        gross_region_vec = fine_to_gross(local_experimentsLocationMatrix, local_humanOntology,local_grossStructures);
        gross_structures_names = local_grossStructures;
    end
    relevant_samples = (gross_region_vec>0);        
    expression = local_experimentsDataMatrix(relevant_samples,:);
    gross_region_vec = gross_region_vec(relevant_samples);
    gene_info = local_gene_info;
    samples2subjects = local_samples2subjects(relevant_samples,:);
    physicalLocation = local_physicalLocation(relevant_samples,:);
end

% ===========================================================
function  [expression, gross_region_vec, gene_info, samples2subjects,gross_structures_names] = load_mouse(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_info;
    persistent local_grossStructures
     
    if isempty(local_expression)
        fprintf('mouse: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'home','lab', 'gal', 'Projects', ...
                                'Limor', 'Data');
        filename = 'ABA_expression_data.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        % Focuse on gross structures
        grossStructures = {'Cerebral cortex'; 'Olfactory areas'; ...
                           'Hippocampal region'; 'Retrohippocampal region'; ...
                           'Striatum'; 'Pallidum'; ...
                           'Hypothalamus'; 'Thalamus'; 'Midbrain'; ...
                           'Cerebellum'; 'Pons'; 'Medulla'};
        num_gross = numel(grossStructures);
        
        inds = cellfind(x.name_of_regions, grossStructures);
        if any(inds==0), error('Name or region not found');end
        expression = x.expression_energy(inds,:);
        gross_region_vec = (1:num_gross);
        
        gene_info.gene_symbols = x.name_of_genes;
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_info = gene_info;
        local_grossStructures = grossStructures;
    end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_info = local_gene_info;
    samples2subjects = ones(size(gross_region_vec)); % the mouse does not have different subjects;
    gross_structures_names = local_grossStructures;
end



% ===========================================================
function  [expression, gross_region_vec,gene_info,samples2subjects,gross_structures_names] = load_kang(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_info;
    persistent local_grossStructures;

    
%     if isempty(local_expression) || ~isempty(parms )
        fprintf('kang: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'cortex','data', 'microarray', 'human', ...
                                'Kang2011');
        filename = 'kang_samples_with_ontology.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        
        % Focuse on gross structures
        %
        % We need to map between 4 layes: 
        % samples->regions->gross->sorted-strurctures
 
        
        % Map from gross to sorted
        grossStructures = { 'Frontal Lobe','Hippocampus',  ...
                            'Occipital Lobe', 'Parietal Lobe', ...
                            'Temporal Lobe', 'Amygdala', 'Striatum', ...
                            'Thalamus', 'Cerebellum'};        
        sorted_inds = cellfind(x.grossRegionNames, grossStructures);
%         [~,sorted_inds] = ismember(grossStructures, x.grossRegionNames);
        assert( all(0 < sorted_inds ) ,'Name of region not found');
        
        % Map from samples to gross: 
        grossRegionSamples = double(x.samples2regions) * double(x.regionOntology');
        samples2grossID = grossRegionSamples * sorted_inds';
        is_sample_good_regions = (samples2grossID>0);
 
      
        % Focuse on adult periods        
        x.periods = x.samples2periods * (1:size(x.samples2periods, 2))';
        
        if isempty(parms)
            ageOfIntrest = {'12-20y', '20-40y', '40-60y', '60+y'};        
            % only adults
        else
            ageOfIntrest = parms;
        end
        
        is_sample_good_ages = ismember(x.period_names(x.periods), ageOfIntrest);
        
        % Focus on adults and on gross regions
        relevent_samples = is_sample_good_regions & is_sample_good_ages';
        expression = x.data(:, relevent_samples)';
        gross_region_vec =  samples2grossID( relevent_samples );
        
        gene_info.gene_symbols = x.gene_names;
        gene_info.entrez_ids = x.entrez;
        
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_info = gene_info;
        local_grossStructures = grossStructures;

%     end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_info = local_gene_info;
    gross_structures_names = local_grossStructures;
    
    samples2subjects = NaN;

end

function  [expression, gross_region_vec,gene_info,samples2subjects,gross_structures_names] = load_kang_all_regions(parms,region_subset)

    fprintf('kang: loading data from disk\n');
    % Load human data 
    data_dirname = fullfile('/', 'cortex','data', 'microarray', 'human', ...
                            'Kang2011');
    filename = 'kang_samples_with_ontology.mat';
    fullname = fullfile(data_dirname, filename);
    x = load(fullname);

    % Get all the sample from regions which are defined in periods 3-15
    
    switch region_subset
        case 'all' 
            sample2region_3_15 = x.samples2regions(:, x.period3_15);           
            gross_structures_names = x.region_names( x.period3_15 );
        case 'cortex'
            cortex_region = {'primary auditory (A1) cortex';'dorsolateral prefrontal cortex';'posterior inferior parietal cortex';'inferior temporal cortex';'primary motor (M1) cortex';'medial prefrontal cortex';'orbital prefrontal cortex';'primary somatosensory (S1) cortex';'superior temporal cortex';'primary visual (V1) cortex';'ventrolateral prefrontal cortex';};
            is_cortex_region = ismember(x.region_names, cortex_region);
            sample2region_3_15 = x.samples2regions(:, x.period3_15 & is_cortex_region);
            gross_structures_names = x.region_names( x.period3_15 & is_cortex_region);
    end
    is_sample_good_regions = any(sample2region_3_15 ,2); 
    sample_region_vector = sample2region_3_15 * (1:size(sample2region_3_15,2))';
    


    % Focuse on adult periods        
    x.periods = x.samples2periods * (1:size(x.samples2periods, 2))';

    if isempty(parms)
        ageOfIntrest = {'12-20y', '20-40y', '40-60y', '60+y'};        
        % only adults
    else
        ageOfIntrest = parms;
    end

    is_sample_good_ages = ismember(x.period_names(x.periods), ageOfIntrest);

    % Focus on adults and on gross regions
    relevent_samples = is_sample_good_regions & is_sample_good_ages';
    expression = x.data(:, relevent_samples)';

    gross_region_vec =  sample_region_vector( relevent_samples );

    gene_info.gene_symbols = x.gene_names;
    gene_info.entrez_ids = x.entrez;

    gross_region_vec = gross_region_vec(:);
    samples2subjects = NaN;

end

% ===========================================================
function  [expression, gross_region_vec,gene_info,samples2subjects,gross_structures_names] = load_kang_cortex(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_info;
    persistent local_grossStructures;

    
    if isempty(local_expression)
        fprintf('kangCortex: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'cortex','data', 'microarray', 'human', ...
                                'Kang2011');
        filename = 'kang_samples_with_ontology.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        % Focuse on gross structures
        %
        % We need to map between 4 layes: 
        % samples->regions->gross->sorted-strurctures
 
        
        % Map from gross to sorted
        grossStructures = { 'Frontal Lobe','Occipital Lobe', 'Parietal Lobe', 'Temporal Lobe'};        
        sorted_inds = cellfind(x.grossRegionNames, grossStructures);
        assert( all(0 < sorted_inds ) ,'Name of region not found');
        
        % Map from samples to gross: 
        truncated_ontology = x.regionOntology(sorted_inds,:);
        grossRegionSamples = double(x.samples2regions) * double(truncated_ontology');
        samples2grossID = grossRegionSamples * ((1:size(grossRegionSamples,2))');
        is_sample_good_regions = (samples2grossID>0);
      
        % Focuse on adult periods        
        x.periods = x.samples2periods * (1:size(x.samples2periods, 2))';
        adults = {'12-20y', '20-40y', '40-60y', '60+y'};        
        is_sample_good_ages = ismember(x.period_names(x.periods), adults);
        
        % Focus on adults and on gross regions
        relevent_samples = is_sample_good_regions & is_sample_good_ages';
        expression = x.data(:, relevent_samples)';
        gross_region_vec =  samples2grossID( relevent_samples );
        
        gene_info.gene_symbols = x.gene_names;
        gene_info.entrez_ids = x.entrez;
        
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_info = gene_info;
        local_grossStructures = grossStructures;

    end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_info = local_gene_info;
    gross_structures_names = local_grossStructures;
    
    samples2subjects = NaN;

end


function  [expression, gross_region_vec,gene_info,samples2subjects,gross_structures_names] = load_brainspan(parms)
% 
    persistent local_expression;
    persistent local_gross_region_vec    
    persistent local_gene_info;
    persistent local_grossStructures
    persistent local_samples2subjects

    if isempty(local_expression)
        fprintf('brainspan: loading data from disk\n');
        % Load human data 
        data_dirname = fullfile('/', 'cortex','data', 'microarray', 'human', ...
                                'brainspan_microarray');
        filename = 'brainspan_microarray.mat';
        fullname = fullfile(data_dirname, filename);
        x = load(fullname);
        
        % Focuse on gross structures
        %
        % We need to map between 4 layes: 
        % samples->regions->gross->sorted-strurctures
 
        
        % Map from gross to sorted
        grossStructures = { 'Frontal Lobe','Hippocampus',  ...
                            'Occipital Lobe', 'Parietal Lobe', ...
                            'Temporal Lobe', 'Amygdala', 'Striatum', ...
                            'Thalamus', 'Cerebellum'};        
        sorted_inds = cellfind(x.grossRegionNames, grossStructures);
        assert( all(0 < sorted_inds ) ,'Name of region not found');
        
        % Map from samples to gross: 
        grossRegionSamples = double(x.samples2regions) * double(x.regionOntology');
        samples2grossID = grossRegionSamples * sorted_inds';
        is_sample_good_regions = (samples2grossID>0);
 
      
        % Focuse on adult periods        
        afterBirth = ~strcmp(x.ages(:,2),'pcw') ;        
        older12 = cell2mat(x.ages(:,1)) >= 12;
        is_sample_good_ages = afterBirth ;%& older12;
        
        % Focus on adults and on gross regions
        relevent_samples = is_sample_good_regions & is_sample_good_ages;
        expression = x.expression(:, relevent_samples)';
        gross_region_vec =  samples2grossID( relevent_samples );
        
        [~,~,samples_subjects] = unique(x.donor_id);
        samples_subjects = samples_subjects(relevent_samples);
        samples2subjects = sparse(1:length(samples_subjects), samples_subjects, ones(size(samples_subjects)) );
        samples2subjects = logical(samples2subjects);
        
        gene_info.gene_symbols = x.gene_symbol;
        gene_info.entrez_ids = x.entrez_id;
        gene_info.ensembl_gene_ids = x.ensembl_gene_id;
        
        local_expression = expression;
        local_gross_region_vec = gross_region_vec(:);
        local_gene_info = gene_info;
        local_grossStructures = grossStructures;
        local_samples2subjects = samples2subjects;

    end
    expression = local_expression;
    gross_region_vec = local_gross_region_vec(:);
    gene_info = local_gene_info;
    gross_structures_names = local_grossStructures;
    samples2subjects = local_samples2subjects;

end


