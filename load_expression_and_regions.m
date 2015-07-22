function [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info, physicalLocation] = load_expression_and_regions(dataname, parms)

physicalLocation = [];

 switch dataname
      case 'human6GrossRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info,physicalLocation] = load_human6(parms);
      case 'human6AllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info,physicalLocation] = load_human6('all regions');
     case 'human6Cortex', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info,physicalLocation] = load_human6_cortex(parms);
     case 'human6CortexAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info,physicalLocation] = load_human6_cortex('all regions');
     case 'kangGrossRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_kang(parms);
     case 'kangAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_kang_all_regions(parms,'all');
    case 'kangCortex', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_kang_cortex(parms);
    case 'kangCortexAllRegions', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_kang_all_regions(parms,'cortex');
    case 'mouse', 
          resolution = parms;
          clear('parms');
          parms.limitRegions = false;
          parms.limitList = {};
          parms.resolution = resolution;
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_mouse(parms);
    case 'mouseLimitRegions', 
        limitList = parms;
        clear('parms');
        parms.limitRegions = true;
        parms.limitList = limitList;
        parms.resolution = 'standard+cortex';
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_mouse(parms);    
    case 'zapalaMouse', 
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_zapala_mouse(parms);
    case 'brainspan',
        [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_info] = load_brainspan(parms);
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
        data_dirname = fullfile('/cortex/data/microarray/human/Hawrylycz_2012');
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

function [expression, gross_region_vec, gene_info, samples2subjects, gross_structures_names] = load_zapala_mouse(parms)

    fprintf('zapala mouse: loading data from disk\n');
    fullname = '/cortex/data/microarray/mouse/Zapala_2005/GSE3594-Zapala-2005.mat';
    x = load(fullname);
    
    % === do some filtering on the samples
    
    %     regions_of_intresets = {'Amygdala';'Bed nucleus of the stria terminalis';'Cerebellum';'Isocortex';'Dentate Gyrus';'Entorhinal Cortex';'Hippocampus';'Hippocampus CA1';'Hippocampus CA3';'Hypothalamus';'Inferior Colliculi';'Medulla';'Thalamus';'Motor Cortex';'Olfactory Bulbs';'Perirhinal Cortex';'Pons';'Striatum';'Superior Colliculi'};
    regions_of_intresets = {'Olfactory Bulbs';'Isocortex';'Motor Cortex';'Entorhinal Cortex';'Perirhinal Cortex';'Hippocampus';'Hippocampus CA1';'Hippocampus CA3';'Dentate Gyrus';'Amygdala';'Bed nucleus of the stria terminalis';'Striatum';'Thalamus';'Hypothalamus';'Inferior Colliculi';'Superior Colliculi'  ;'Cerebellum';'Pons';'Medulla'};
    
    regions_not_included = {'Adrenal';'Brown Adipose Tissue';'Choroid Plexus';'Heart';'Kidney';'Liver';'Muscle';'PAG';'Pituitary';'Retina';'Spinal Cord';'Spleen';'Testes';'Thymus';'White Adipose Tissue'};
    % === done with the filtering
    
    expression = x.gene_expression_data';

    gene_info.gene_symbols = x.gene_symbol;
    gene_info.entrez_ids = x.entrez;
    gene_info.gene_full_name = x.gene_title;
    gene_info.probe_id = x.gpl_probe_id;
        
    gross_structures_names = regions_of_intresets;
    
    % Change the name of BNST to 'Bed nucleus of the stria terminalis'
    sample_sources = x.sample_sources;
    sample_sources(ismember(sample_sources,'BNST') ) = {'Bed nucleus of the stria terminalis'};
    sample_sources(ismember(sample_sources,'Midbrain-Thalamus') ) = {'Thalamus'};
    sample_sources(ismember(sample_sources,'Midbrain-thalamus') ) = {'Thalamus'};
    sample_sources(ismember(sample_sources,'Cortex') ) = {'Isocortex'};

    
    [~ ,gross_region_vec] = ismember(sample_sources, gross_structures_names);
    non_relevent_samples = gross_region_vec == 0;
    gross_region_vec( non_relevent_samples ) = [];
    expression(non_relevent_samples,:) = [];
    
    samples2subjects = ones(size(gross_region_vec)); % the mouse does not have different subjects;
    
    [gross_region_vec, resort_ind] = sort(gross_region_vec);
    expression = expression(resort_ind,:);
end
% ===========================================================
% function  [expression, gross_region_vec, gene_info, samples2subjects,gross_structures_names] = load_mouse(parms)
% % 
%     persistent local_expression;
%     persistent local_gross_region_vec    
%     persistent local_gene_info;
%     persistent local_grossStructures
%      
%     if isempty(local_expression)
%         fprintf('mouse: loading data from disk\n');
%         % Load human data 
%         data_dirname = fullfile('/', 'home','lab', 'gal', 'Projects', ...
%                                 'Limor', 'Data');
%         filename = 'ABA_expression_data.mat';
%         fullname = fullfile(data_dirname, filename);
%         x = load(fullname);
%         
%         % Focuse on gross structures
%         grossStructures = {'Cerebral cortex'; 'Olfactory areas'; ...
%                            'Hippocampal region'; 'Retrohippocampal region'; ...
%                            'Striatum'; 'Pallidum'; ...
%                            'Hypothalamus'; 'Thalamus'; 'Midbrain'; ...
%                            'Cerebellum'; 'Pons'; 'Medulla'};
%         num_gross = numel(grossStructures);
%         
%         inds = cellfind(x.name_of_regions, grossStructures);
%         if any(inds==0), error('Name or region not found');end
%         expression = x.expression_energy(inds,:);
%         gross_region_vec = (1:num_gross);
%         
%         gene_info.gene_symbols = x.name_of_genes;
%         local_expression = expression;
%         local_gross_region_vec = gross_region_vec(:);
%         local_gene_info = gene_info;
%         local_grossStructures = grossStructures;
%     end
%     expression = local_expression;
%     gross_region_vec = local_gross_region_vec(:);
%     gene_info = local_gene_info;
%     samples2subjects = ones(size(gross_region_vec)); % the mouse does not have different subjects;
%     gross_structures_names = local_grossStructures;
% end



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
    
    [~,~,subjectId] = unique(x.subject_id);
    samples2subjects = logical(full(sparse(1:length(subjectId) , subjectId, ones(size(subjectId)) ) ));
    samples2subjects = samples2subjects(relevent_samples,:);
    samples2subjects = samples2subjects(:, any(samples2subjects,1) ); %remove all subjects which does not have any samples representing them.
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
    
    [~,~,subjectId] = unique(x.subject_id);
    samples2subjects = logical(full(sparse(1:length(subjectId) , subjectId, ones(size(subjectId)) ) ));
    samples2subjects = samples2subjects(relevent_samples,:);
    samples2subjects = samples2subjects(:, any(samples2subjects,1) ); %remove all subjects which does not have any samples representing them.

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
    
    [~,~,subjectId] = unique(x.subject_id);
    samples2subjects = logical(full(sparse(1:length(subjectId) , subjectId, ones(size(subjectId)) ) ));
    samples2subjects = samples2subjects(relevent_samples,:);
    samples2subjects = samples2subjects(:, any(samples2subjects,1) ); %remove all subjects which does not have any samples representing them.

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


function [expression, gross_region_vec, gene_info,samples2subjects,gross_structures_info] = load_mouse(parms)
    %+=====================================================
    % region_resolution can be:  'standard'    'cortex'    'standard+cortex'    'fine'    'big12'    'cortexLayers'
    %
    % should return:
    % expression    samples X gene matrix
    % gross_region_vec      samples X 1  , with indices for the region (1...N)
    % gene_info -  a struct with fields entrez_id, gene_symbols, gene_ids
    % samples2subjects - not relevent should be a one vector
    % gross_structures_names - should be a cell array with columns:    id, atlas_id, ack, region_name
    %======================================================

    region_resolution = parms.resolution;

    fprintf('Allen mouse: loading data from disk\n');
    addpath('/home/lab/lior/Projects/buildStructureOntology/');
    load('/home/lab/lior/Projects/buildStructureOntology/mouseOntologyObject.mat')

    addpath(genpath('/cortex/data/ISH/mouse/AllenBrainToolbox' ));
    load( '/cortex/data/ISH/mouse/AllenBrainToolbox/atlasData/refAtlas.mat' );
    dbstop if error;
    load( '/cortex/data/ISH/mouse/AllenBrainToolbox/atlasData/ExpEnergy.mat' );
    cor = Ref.Coronal;

    gene_info.gene_symbols = get_genes( cor, 'allNoDup', 'allen' );
    gene_info.entrez_ids = get_genes( cor, 'allNoDup', 'entrez' );

    ann = cor.Annotations;
    resolution_index = find( strcmp( region_resolution, ann.identifier) );
    
    annotions3D =  get_annotation ( cor , region_resolution );
    gross_structures_ids = mat2cell( ann.ids{ resolution_index } ,ones(1,size( ann.ids{ resolution_index } ,1)),1);
%     gross_structures_ids = cellfun(@num2str, gross_structures_ids, 'UniformOutput',false);

    mouse_3d_symbols = ann.symbols{ resolution_index } ;
    mouse_3d_names = ann.labels{ resolution_index } ;
    gross_region_vec_orig_indices = ann.classification{ resolution_index };
    expression = getSampleInFilter(D, resolution_index, cor, ann);
    region_indices = ann.ids{ resolution_index };
    

    [is_ontology_member, ind_in_ontolgy] =    ismember( mouse_3d_symbols , mouseOntology.structureLabels(:,3) ) ;
    if( ~all(is_ontology_member) )
        fprintf('not all regions from the 3D mouse data are in the mouse ontlogy\n');
        disp(ann.labels{resolution_index}(~is_ontology_member));

        % limit the sample and region to only those in the ontologyremove...
        ind_in_ontolgy( ind_in_ontolgy ==0 ) = [];
        region_indices(~is_ontology_member) = [];
        mouse_3d_symbols(~is_ontology_member) = [];
        mouse_3d_names(~is_ontology_member) = [];
    end
    
    gross_structures_info = mouseOntology.structureLabels(ind_in_ontolgy, :);
    [is_sample_a_member, gross_region_vec] = ismember(gross_region_vec_orig_indices,  region_indices  );
    if (~all(is_sample_a_member))
        fprintf('some samples does not have an associated region\n');
        gross_region_vec(~is_sample_a_member) = [];
        expression(~is_sample_a_member,:) = [];
    end

    samples2subjects = ones(size(gross_region_vec));

    if parms.limitRegions
        
        [gross_region_vec, gross_structures_info] = mapFineToGross(gross_region_vec, parms.limitList , gross_structures_info(:,4), mouseOntology);
        
        no_parent_regions = gross_region_vec ==0;
        expression = expression(~no_parent_regions, :);
        gross_region_vec = gross_region_vec(~no_parent_regions);
        samples2subjects = samples2subjects(~no_parent_regions);
        
        regions_with_samples = unique(gross_region_vec);
        region_has_samples = false(size(gross_structures_info,1),1 ); region_has_samples(regions_with_samples) = true;
        if ~all(region_has_samples)
            fprintf('The following regions do not have samples associates with it:\n');
            disp( gross_structures_info(~region_has_samples,4) );
        end
        
    end
end

function [gross_region_vec, coarseRegionsInfo] = mapFineToGross(regionPerSample, coarseRegions, fineRegions, ontology)
% if inclusive is true, every fine region region which is a part of a
% coarse region is deleted, every fine region which is not a part of a
% coarse region is left as is.
% if inclusive is false, fine regions are mapped to coarse regions,
% fine regions which does not have a mapping are ommited.


    % create a belongs to matrix with dim  (num_samples X num_fine_regions) 
    sample_region_matrix = full(sparse(1:length(regionPerSample), regionPerSample, ones(size(regionPerSample)) )) ;
    [num_samples, num_fine_regions] = size(sample_region_matrix);
    num_regions_in_ontology = size(ontology.structureLabels,1);
    assert(num_fine_regions == length(fineRegions),  'the maximum index in regionPerSample should be of length(fineRegions)');

    [is_member, fine_indices] = ismember(fineRegions, ontology.structureLabels(:,4) );
    assert(all(is_member), 'there is a region in the fine-regions list which is not in the onotlogy');

    [is_member, coarse_indices] = ismember(coarseRegions, ontology.structureLabels(:,4) );
    assert(all(is_member), 'there is a region in the coarse-regions list which is not in the onotlogy');

    assert(all(strcmp( fineRegions , ontology.structureLabels(fine_indices,4)  )), 'there is a region in the fine-regions list which is not in the onotlogy');

    isChildMatrix = ontology.allChildNodes();
    coarseRegionsInfo = ontology.structureLabels(coarse_indices,:);
    isChildMatrix = isChildMatrix(coarse_indices, fine_indices) ;

    sample_region_matrix_gross = sample_region_matrix * isChildMatrix';
    if( ~all( sum(isChildMatrix,1) == ones(1, size(isChildMatrix,2) ) ))
        fprintf('The regions list is not comprehensive. Some regions and samples are not in the new list\n');
    end
    assert( max( sum(sample_region_matrix_gross,2)) <= 1 ,'each sample should have at most one region');


    one_to_n = (1:size(sample_region_matrix_gross,2) )';
    gross_region_vec = sample_region_matrix_gross * one_to_n;
end

function DFiltered = getSampleInFilter(D, identifierIndex, cor, ann)


    % Compute the s e t o f rows us ing a volume o f
    %indices and the filter2.2. FROM MATRICES TO VOLUMES 11

    brainFilter = get_voxel_filter ( cor , 'brainVox' );
    numVox = numel ( brainFilter );
    % l a b e l the voxe ls by in t e g e r s
    indsBrainVoxels = 1 : numVox ;
    % arrange the in t e g e r s in a volume
     indsVol = make_volume_from_labels ( indsBrainVoxels , brainFilter );
     filter = get_voxel_filter ( cor , ann.filter{ identifierIndex } );
     % r e s t r i c t the volume to the voxe ls that are in the f i l t e r
     indsFiltered = indsVol ( filter );
     DFiltered = D( indsFiltered , : );

end
