function empirical_pvalue = getEmpiricalPvalues(scores, randomScores)
    randomScores = randomScores(:);
    scores = scores(:);

    
%%=============   Method I (Slow) ===============
% %for each value in scores loop over all values in random score and count
% % how many random values are bigger than the score.
%
%     empirical_pvalue = nan(size(scores));
%     for i =1:length(scores)
%         empirical_pvalue(i) = sum(  randomScores > scores(i) );
%     end
%     empirical_pvalue = (empirical_pvalue +1)/ (length(randomScores) +1);


%%=============   Method II (Faster) ===============
% Do the same just use one line with arrayfun.
% % how many random values are bigger than the score.
%
% empirical_pvalue = arrayfun(@(x) sum( randomScores > x ), scores);
% empirical_pvalue = (empirical_pvalue +1)/ (length(randomScores) +1);


%%=============   Method III (Fastest) ===============
% Find the rank of each value in scores inside the joint vector - [scores,random]
% Than reduce from this ranking the ranking of the vector scores
    together = [scores ; randomScores];
    together_ordered = tiedrank(together);
    together_ordered_just_scores = together_ordered(1:length(scores));

    scores_ordered = tiedrank(scores);
    empirical_pvalue = together_ordered_just_scores  - scores_ordered;
    empirical_pvalue = length(randomScores) - empirical_pvalue;

    empirical_pvalue = (empirical_pvalue +1)/ (length(randomScores) +1);

end

function geneLargerThanX = moreThenXPercent(scores, randomScores, x_percent)

randomScoreSorted = sort(randomScores(:) );

index = round( x_percent * length(randomScoreSorted) );
geneLargerThanX = scores > randomScoreSorted(index);
how_many_larger = sum(  geneLargerThanX );

fprintf( '%d%% (%d) of the genes are larger \n', round( how_many_larger / length(scores)*100), how_many_larger);
end