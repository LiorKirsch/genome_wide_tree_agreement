function drawBars(values,std,labels,colors,yAxisLabel)
% draw bars with a differnt color and label for each bar
    set(gca, 'FontSize',16);
    hold on;
    
    for i=1:length(values)
        bar(i,values(i),'FaceColor',colors(i,:));
        errorbar(values,std,'.');
    end
    
    hold off;
    
    ax = gca;
    ax.XTick = 1:length(labels);
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 45;
%     xticklabel_rotate(1:length(values),40,labels, 'fontSize',14) ;
    ylabel(yAxisLabel);
%     figureHandle = gcf;
%     set(findall(figureHandle,'type','text'),'fontSize',14);
    
end
