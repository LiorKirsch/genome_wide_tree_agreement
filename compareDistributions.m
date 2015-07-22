function compareDistributions(distA, distB, nameA, nameB)
    pvalue = ranksum(distA, distB);
    meadianA = median(distA);
    meadianB = median(distB);
    fprintf('%s (median %g) - %s (median %g):\t\t %g both', nameA,meadianA,nameB,meadianB, pvalue);
    pvalueright = ranksum(distA, distB,'tail','right');
    fprintf(',\t\t %g right', pvalueright);
    pvalueleft = ranksum(distA, distB,'tail','left');
    fprintf(',\t\t %g left', pvalueleft);
    fprintf('\n');
 
end