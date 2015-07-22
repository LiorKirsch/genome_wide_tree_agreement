function draw_groups_dist(name,cell_scores, legend_strings, colors, xlimit,ylimit, open_new_fig)
if ~exist('open_new_fig','var')
    open_new_fig = true;
end

    %======= some defaults ====
    line_width = 3;
    x_label =  'BRO-agreement (ABA6-2013)';
    %==========================
    switch name
        case 'Human6' 
            figure_header = 'Human6';
            file_name = 'dist_random';
            do_ylim= false;
            yticks =[0:0.2:1];
        case 'Kang'
            figure_header = 'Kang';
            file_name = 'dist_kang_random';
            x_label = 'BRO-agreement (Kang-2011)';
            do_ylim= false;
            yticks =[0:0.2:1];
        case 'Human6 cahoy'
            figure_header = 'Human6 cahoy';
            file_name = 'dist_cellType';
            do_ylim= true;
        case 'Human6 housekeeping'
            figure_header = 'Human6 housekeeping';
            file_name = 'dist_house';
            do_ylim= true;
        case 'Human6 axon guidance'
            figure_header = 'Human6 housekeeping';
            file_name = 'dist_axon';
            do_ylim= true;
        case 'Human6 serotonin'
            figure_header = 'Human6 serotonin';
            file_name = 'dist_serotonin';
            do_ylim= true;
        case 'Age'
            figure_header = 'Age';
            file_name = 'dist_age';
            do_ylim= true;
        otherwise
            error('unsupported figure type - %s',name);
                   
    end
    
%     for i = 1:length(legend_strings)
%         legend_strings{i} = sprintf('%s (%d)', legend_strings{i}, length(cell_scores{i}) );
%     end

if open_new_fig
    figure('name',sprintf('%s distribution', figure_header));  
end

ploth = displayMultiDist(cell_scores, x_label,legend_strings,xlimit);
set(ploth,'LineWidth',  line_width );
for i =1:size(colors,1)
    set(ploth(i),'Color', colors(i,:) ); 
end

redoDistImage(xlimit); set(ploth,'LineWidth',  line_width );

if do_ylim
    ylim(ylimit);
end
if exist('yticks', 'var')
    set(gca,'ytick',yticks);
end
saveFigure(gcf, fullfile('figures',file_name), 'png');
saveFigure(gcf, fullfile('figures',file_name), 'tiff');
saveFigure(gcf, fullfile('figures',file_name), 'eps');


end

function ploth = displayMultiDist(scores_cell, xlabel_string,legendStrings,range)
    spacing = linspace(range(1),range(2),50);
    vals_in_bins = nan(length(spacing), length(scores_cell));
    for i =1 :length(scores_cell);
        scores = scores_cell{i};
        vals_in_bins(:,i) = histc(scores(:), spacing);
%         vals_in_bins(:,i) = smooth(vals_in_bins(:,i), 0.1, 'moving');
        vals_in_bins(:,i) = vals_in_bins(:,i) / sum(vals_in_bins(:,i));
    end

    vals_in_bins_smooth = smooth(vals_in_bins,0.02);
    vals_in_bins_smooth = reshape(vals_in_bins_smooth,size(vals_in_bins));
    
    ploth = plot(spacing,  vals_in_bins_smooth    );
    xlabel(xlabel_string,'fontsize',20);
    h_legend = legend(legendStrings);
    legend('boxoff');
    set(h_legend,'FontSize',20);
    set(gca,'box','off');
    set(gca,'Fontsize',17);
    set(ploth,'LineWidth',  2 );
end

function redoDistImage(xlimit)
    xtick = -0.2:0.2:1;
    ytick = -0.2:0.1:1;

    set(gca,'Fontsize',20);         %set(gca,'xscale','log');
    set(gca,'ytick',ytick); set(gca,'xtick',xtick);
    xlim( xlimit)
end
