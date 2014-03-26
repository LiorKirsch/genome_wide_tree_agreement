function empirical_pvalue = getEmpiricalPvalues(scores, randomScores)
    randomScores = randomScores(:);
    scores = scores(:);

    empirical_pvalue = arrayfun(@(x) sum( randomScores > x ), scores);
    
%     empirical_pvalue = nan(size(scores));
%     for i =1:length(scores)
%         empirical_pvalue(i) = sum(  randomScores > scores(i) );
%     end

    empirical_pvalue = empirical_pvalue / length(randomScores);

end