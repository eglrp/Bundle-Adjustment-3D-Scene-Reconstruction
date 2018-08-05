function [newS, newCamMat]=rebuild_camera_matrix(S, newCamMat, size_c, opts)
    if size(newCamMat, 1)>0
        G=newCamMat;
    else
        G=build_adjacency_matrix(S, size_c);
        if opts.reducedCameraMatrixType==0
            % Do nothing
        elseif opts.reducedCameraMatrixType==1
            G=diag(ones(size(G, 1), 1));
        elseif opts.reducedCameraMatrixType==2
            [G, discovered_vertices]=build_max_leaf_spanning_tree(G);        
            G=G+diag(ones(size(G, 1), 1));
        elseif opts.reducedCameraMatrixType==3
            [T, discovered_vertices]=build_max_leaf_spanning_tree(G);
            G=build_skeletal_graph(G, T, discovered_vertices, opts.spanner_t);             
            G=G+diag(ones(size(G, 1), 1));
        end
        newCamMat=G;
    end
    
    [row, col]=find(G);
    newS=spalloc(size(S, 1), size(S, 2), size(row, 1)*size_c*size_c);
    for i=1:size(row, 1)
        rowBlock=(row(i)-1)*size_c+1:row(i)*size_c;
        colBlock=(col(i)-1)*size_c+1:col(i)*size_c;
        newS(rowBlock, colBlock)=S(rowBlock, colBlock); %#ok<SPRIX>
    end    
end

function G=build_adjacency_matrix(S, size_c)
    [row, col]=find(S);
    row=ceil(row/size_c);
    col=ceil(col/size_c);
    G=sparse(row, col, ones(size(row, 1), 1));
    G=G-diag(diag(G));
end

function G=process_graph(G)
    [T, discovered_vertices]=build_max_leaf_spanning_tree(G);    
    spanner_t=2;
    G=build_skeletal_graph(G, T, discovered_vertices, spanner_t);    
    G=G+diag(ones(size(G, 1), 1));
end

function [T, discovered_vertices]=build_max_leaf_spanning_tree(G)      
    root_vtx=-1;
    max_degree=-1;
    for i=1:size(G, 1)
        degree=size(find(G(i, :)), 2);
        if max_degree<degree
            max_degree=degree;
            root_vtx=i;
        end        
    end    
    saved_edges=[];
    discovered_vertices=zeros(size(G, 1), 1);  
    discovered_vertices(root_vtx)=1;
    process_vertex(root_vtx);    
    % Builds the max-leaf spanning tree, also deletes these edges from G
    function process_vertex(vtx)
        nbrs=find(G(vtx, :));
        cur_vertices=[];        
        for k=1:size(nbrs, 2)
            nbr=nbrs(k);
            if discovered_vertices(nbr)==0
                saved_edges=[saved_edges; vtx nbr]; %#ok<AGROW>                                
                cur_vertices=[cur_vertices; nbr]; %#ok<AGROW>
                discovered_vertices(nbr)=1;
            end
        end
        children_max_degree=-1;
        cand_vtx=-1;
        for k=1:size(cur_vertices, 1)
            next_vtx=cur_vertices(k);
            all_nbrs=find(G(next_vtx, :));
            degree=size(find(discovered_vertices(all_nbrs)==0), 1);
            if children_max_degree<degree
                children_max_degree=degree;
                cand_vtx=next_vtx;
            end
        end
        % vtx is an interior node, not a tree-leaf
        if cand_vtx>0
            process_vertex(cand_vtx);
            discovered_vertices(vtx)=2;
        end       
    end
    T=sparse(size(G, 1), size(G, 2));
    for i=1:size(saved_edges, 1)        
        T(saved_edges(i, 1), saved_edges(i, 2))=1;
        T(saved_edges(i, 2), saved_edges(i, 1))=1;
    end
end

function T=build_skeletal_graph(G, T, discovered_vertices, spanner_t)
    interior_nodes=find(discovered_vertices==2);
    leaves=find(discovered_vertices==1);
    edge_list=[];
    for i=1:size(interior_nodes, 1)
        for j=i+1:size(interior_nodes, 1)
           if(G(interior_nodes(i), interior_nodes(j))~=0)&&(T(interior_nodes(i), interior_nodes(j))==0)
               edge_list=[edge_list; interior_nodes(i) interior_nodes(j)];
           end
        end
        
        for j=1:size(leaves, 1)
            if(G(interior_nodes(i), leaves(j))~=0)&&(T(interior_nodes(i), leaves(j))==0)
                edge_list=[edge_list; interior_nodes(i) leaves(j)];
            end
        end
    end
    for i=1:size(leaves, 1)
        for j=i+1:size(leaves, 1)
            if(G(leaves(i), leaves(j))~=0)&&(T(leaves(i), leaves(j))==0)
                edge_list=[edge_list; leaves(i) leaves(j)];
            end
        end
    end
    for i=1:size(edge_list, 1)        
        pathLength=graphshortestpath(T, edge_list(i, 1), edge_list(i, 2), 'Directed', false, 'Method', 'BFS');
        if pathLength>spanner_t
            T(edge_list(i, 1), edge_list(i, 2))=1;
            T(edge_list(i, 2), edge_list(i, 1))=1;
        end
    end    
end

