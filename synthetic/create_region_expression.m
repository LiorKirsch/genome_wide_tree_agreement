function sample_expression  = create_region_expression( tree_node_expression, region_vector, conf)
   
%=============== using for loop =============
%     sample_expression = nan(size(region_vector,1),1);
%     for i = 1:length(region_vector)
%         current_region_index = region_vector(i);
%         
%         random_noise = rand(1)*2 -1;
%         random_noise = random_noise * conf.sample_noise;
% 
%         sample_expression(i) = tree_node_expression(current_region_index) + random_noise;
%     end
%============================================
    expressionDim = size(tree_node_expression,2);

    random_noise = rand(length(region_vector),expressionDim)*2 -1;
    random_noise = random_noise * conf.sample_noise;
    sample_expression = tree_node_expression(region_vector,:) + random_noise;
end