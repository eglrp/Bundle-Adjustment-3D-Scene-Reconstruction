function isForest=checkIsForest(T)    
    n=size(T, 1);
    discovered_vertices=zeros(n, 1);
    unassigned_vertices=find(discovered_vertices==0);
    isForest=1;
    while size(unassigned_vertices, 1)>0
        root=unassigned_vertices(1);
        dists=graphshortestpath(T, root, 'Directed', false, 'Method', 'BFS');
        treeNodes=[];
        for i=1:n
            if dists(i)~=Inf
                discovered_vertices(i)=1;
                treeNodes=[treeNodes; i];
            end
        end
        T2=T(treeNodes, treeNodes);
        if nnz(T2)~=2*(size(treeNodes, 1)-1)
            disp('bad tree');
            isForest=0;
        end
        unassigned_vertices=find(discovered_vertices==0);
    end    
end