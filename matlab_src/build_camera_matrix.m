function G=build_camera_matrix(G, opts)    
    if opts.reducedCameraMatrixType==0
        G=G+diag(ones(size(G, 1), 1));        
    elseif opts.reducedCameraMatrixType==1
        G=diag(ones(size(G, 1), 1));
    elseif opts.reducedCameraMatrixType==2
        [G, ~]=build_max_leaf_spanning_tree(G); 
        if ~checkIsForest(G)
            disp('T is not a forest');
        end
        G=G+diag(ones(size(G, 1), 1));
    elseif opts.reducedCameraMatrixType==3
        [T, discovered_vertices]=build_max_leaf_spanning_tree(G);
        disp('Finish MLST construction');
        if ~checkIsForest(T)
            disp('T is not a forest');
        end
        G=build_skeletal_graph(G, T, discovered_vertices, opts.spanner_t);             
        G=G+diag(ones(size(G, 1), 1));
    end     
end

function [T, discovered_vertices]=build_max_leaf_spanning_tree(G)
    discovered_vertices=zeros(size(G, 1), 1);
    unassigned_vertices=find(discovered_vertices==0);   
    saved_edges=[];
    while size(unassigned_vertices, 1)>0
        [saved_edges, discovered_vertices]=build_single_max_leaf_spanning_tree(G, unassigned_vertices, discovered_vertices, saved_edges);
        unassigned_vertices=find(discovered_vertices==0);
    end
    T=sparse(size(G, 1), size(G, 2), size(G, 1)*2);
    for i=1:size(saved_edges, 1)        
        T(saved_edges(i, 1), saved_edges(i, 2))=1;
        T(saved_edges(i, 2), saved_edges(i, 1))=1;
    end
    T=T-diag(diag(T));
end

function [saved_edges, discovered_vertices]=build_single_max_leaf_spanning_tree(G, unassigned_vertices, discovered_vertices, saved_edges)      
    root_vtx=-1;
    max_degree=-1;
    for i=1:size(unassigned_vertices, 1)
        vtx=unassigned_vertices(i);
        degree=size(find(G(vtx, :)), 2);
        if max_degree<degree
            max_degree=degree;
            root_vtx=vtx;
        end        
    end    
    
    cand_vertices=root_vtx;
    unmarked_deg=zeros(size(G, 1), 1);
    for i=1:size(G, 1)
        unmarked_deg(i)=size(find(G(i, :)), 2);
    end
    while numel(cand_vertices)>0
        max_unmarked_deg=-1;
        next_vtx=-1;
        for i=1:size(cand_vertices, 1)            
            vtx=cand_vertices(i);            
            if unmarked_deg(vtx)>max_unmarked_deg
                max_unmarked_deg=unmarked_deg(vtx);
                next_vtx=i;
            end                        
        end
        expand_vertex(cand_vertices(next_vtx));
        cand_vertices(i)=[];
    end
    
    function expand_vertex(vid)
        discovered_vertices(vid)=1;
        nbrs=find(G(vid, :));
        for k=1:size(nbrs, 2)
            nbr_idx=nbrs(k);
            if discovered_vertices(nbr_idx)==0
                unmarked_deg(nbr_idx)=unmarked_deg(nbr_idx)-1;
            end
        end
        for k=1:size(nbrs, 2)
            nbr=nbrs(k);
            if discovered_vertices(nbr)==0
                saved_edges=[saved_edges; vid nbr];
                discovered_vertices(nbr)=1;
                cand_vertices=[cand_vertices; nbr];
                nbr_nbrs=find(G(nbr, :));
                for l=1:size(nbr_nbrs, 2)
                    nbr_nbr_idx=nbr_nbrs(l);
                    if discovered_vertices(nbr_nbr_idx)==0
                        unmarked_deg(nbr_nbr_idx)=unmarked_deg(nbr_nbr_idx)-1;
                    end
                end
            end
        end
        discovered_vertices(vid)=2;
    end
    
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
            degree=sum(discovered_vertices(all_nbrs)==0);
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
end

function T=build_skeletal_graph(G, T, discovered_vertices, spanner_t)
    interior_nodes=find(discovered_vertices==2);
    leaves=find(discovered_vertices==1);
    m=size(interior_nodes, 1);
    n=size(leaves, 1);
    edge_list=zeros(m*(m-1)/2+m*n+n*(n-1)/2, 2);
    edge_list_counter=1;
    for i=1:m
        for j=i+1:m
           if(G(interior_nodes(i), interior_nodes(j))~=0)&&(T(interior_nodes(i), interior_nodes(j))==0)
               edge_list(edge_list_counter, :)=[interior_nodes(i) interior_nodes(j)];
               edge_list_counter=edge_list_counter+1;
           end
        end
        
        for j=1:n
            if(G(interior_nodes(i), leaves(j))~=0)&&(T(interior_nodes(i), leaves(j))==0)
                edge_list(edge_list_counter, :)=[interior_nodes(i) leaves(j)];
                edge_list_counter=edge_list_counter+1;
            end
        end
    end
    for i=1:n
        for j=i+1:n
            if(G(leaves(i), leaves(j))~=0)&&(T(leaves(i), leaves(j))==0)
                edge_list(edge_list_counter, :)=[leaves(i) leaves(j)];
                edge_list_counter=edge_list_counter+1;
            end
        end
    end
    for i=1:size(edge_list, 1)
        if edge_list(i, 1)==0 && edge_list(i, 2)==0
            break;
        end
        pathLength=graphshortestpath(T, edge_list(i, 1), edge_list(i, 2), 'Directed', false, 'Method', 'BFS');
        if pathLength>spanner_t
            T(edge_list(i, 1), edge_list(i, 2))=1;
            T(edge_list(i, 2), edge_list(i, 1))=1;
        end
    end    
end