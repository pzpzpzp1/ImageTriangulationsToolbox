% Draws N x 1 elements based on edge score.
% First draws treating score like a distribution. (with replacement)
% Then draws to make up for previously redundant draws using greedy sorting of the score. 
% modifcation to score ensures the second stage will not draw anything already drawn from the first stage.
% Might return less than N elements but ONLY if there are actually fewer edges than N.
% returns edge indices ordered by score
function edgeInds = drawEdgesToSplit(N, score)
    nE = numel(score);
    if N >= nE
        edgeInds = [1:nE]';
        return;
    end
    
    % Draw according to pdf with replacement. 
    % This should help not let the edges to divide get too clustered.
    % edgeInds = unique(randsample(1:nE, N, true, score));
    edgeInds = unique(randsample(1:nE, floor(N/10), true, score));
    
    % draw the remainder from N greedily
    score2 = score + eps;
    score2(edgeInds)=0;
    [~,perm] = sort(score2,'desc');
    edgelist = 1:nE;
    priorityEdgeList = edgelist(perm);
    remainingN = N - numel(edgeInds);
    edgeInds = unique([priorityEdgeList(1:remainingN) edgeInds])';
    
    % return edgeInds sorted by score
    [~,perm] = sort(score(edgeInds),'desc');
    edgeInds = edgeInds(perm);
end