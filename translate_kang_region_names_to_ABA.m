function [region_name_aba, region_sym_aba] = translate_kang_region_names_to_ABA(region_names)

kang_region_names = {
 'primary auditory (A1) cortex';
 'amygdala';
 'cerebellar cortex';
 'dorsolateral prefrontal cortex';
 'hippocampus';
 'posterior inferior parietal cortex';
 'inferior temporal cortex';
 'primary motor (M1) cortex';
 'mediodorsal nucleus of the thalamus';
 'medial prefrontal cortex';
 'orbital prefrontal cortex';
 'primary somatosensory (S1) cortex';
 'superior temporal cortex';
 'striatum';
 'primary visual (V1) cortex';
 'ventrolateral prefrontal cortex';
 
  'Frontal Lobe';
  'Hippocampus';
  'Occipital Lobe';
  'Parietal Lobe';
  'Temporal Lobe';
  'Thalamus';
  'Cerebellum';
 };


kang_region_translation = {
'A1C','primary auditory cortex (core)';
'AMY','amygdaloid complex';
'CBC','cerebellar cortex';
'DFC','dorsolateral prefrontal cortex';
'HIP','hippocampus (hippocampal formation)';
'IPC','posteroventral (inferior) parietal cortex';
'ITC','inferolateral temporal cortex (area TEv, area 20)';
'M1C','primary motor cortex (area M1, area 4)';
'MD','mediodorsal nucleus of thalamus';
'A10','frontal polar cortex (area 10)';
'OFC','orbital frontal cortex';
'S1C','primary somatosensory cortex (area S1, areas 3,1,2)';
'SLTC','superolateral temporal cortex';
'STR','striatum';
'V1C','primary visual cortex (striate cortex, area V1/17)';
'VFC','ventrolateral prefrontal cortex';


 'FroL','Frontal Lobe'; % 'frontal lobe';
 'HIP','hippocampal formation';  % 'hippocampus (hippocampal formation)';
 'OccL','Occipital Lobe'; % 'occipital lobe';
 'ParL','Parietal Lobe'; % 'parietal lobe';
 'TemL','Temporal Lobe'; % 'temporal lobe';
 'THM','Thalamus'; % 'thalamus';
 'CB','Cerebellum'; % 'cerebellum'; 
};

    region_name_aba = cell(length(region_names),1);
    region_sym_aba = cell(length(region_names),1);
    for i =1:length(region_names)
       region_index = strcmp(region_names{i} ,  kang_region_names) ;
       assert( sum(region_index) ==1, sprintf('each region name should have exact 1 translation (error in region %d)', i));
       translation = kang_region_translation(region_index,:);
       region_sym_aba{i} = translation{1};
       region_name_aba{i} = translation{2};
    end
end
