function [distsTable, dists]=debug(G)
    [Grow, Gcol]=find(G);
    G=sparse(Grow, Gcol, ones(size(Grow)));
    distsTable=zeros(size(G));
    for i=1:size(G, 1)
        dists=graphshortestpath(G, i);
        for j=1:size(G, 1)
            distsTable(i, j)=dists(j);
        end
    end
    n=size(G, 1);
    dists=zeros(n*(n-1)/2, 1);
    ctr=1;
    for i=1:size(G, 1)
        for j=1:i-1
            dists(ctr)=distsTable(i, j);
            ctr=ctr+1;
        end
    end
end

