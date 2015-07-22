function printGeneNamesOrdered(geneScore, geneSym,fileName,direction)
% direction should be 'descend' or 'ascend'

    [sortedScore, sortInd] = sort(geneScore, direction);
    fid = fopen(fileName, 'w');
    sortedGeneSym = geneSym(sortInd);
    for i =1:length(sortedGeneSym)
        fprintf(fid,'%g,%s\r\n', sortedScore(i), sortedGeneSym{i} );
    end
    fclose(fid);
end