% f is the objective function value with input state_x vector
% Outputs the Jacobian in the form J=[B E; E' C] and the gradient vector
% g=J'*f
function [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D,opts)    
    if isfield(opts, 'obsMaskPrecon') && isfield(opts, 'obsMask')
        [B, C, E, B2, C2, E2, g, g2, Bindices, Cindices, Eindices, B2indices, C2indices, E2indices]=calculate_jacobian_croute(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.n_pts, bundleData.cidx, bundleData.pidx, D, opts.obsMaskPrecon, size(opts.obsMaskPrecon, 1), opts.obsMask, size(opts.obsMask, 1), f, f2);
        B2=sparse(B2indices(1:2:size(B2indices, 1)-1), B2indices(2:2:size(B2indices, 1)), B2);
        clear B2indices;
        C2=sparse(C2indices(1:2:size(C2indices, 1)-1), C2indices(2:2:size(C2indices, 1)), C2);
        clear C2indices;
        invalidE2entries=find(E2indices==0);
        E2indices(invalidE2entries)=[];
        E2(invalidE2entries(2:2:size(invalidE2entries, 1))/2)=[];
        E2=sparse(E2indices(1:2:size(E2indices, 1)-1), E2indices(2:2:size(E2indices, 1)), E2);
        clear E2indices;
        diag_weights=0;
        B=sparse(Bindices(1:2:size(Bindices, 1)-1), Bindices(2:2:size(Bindices, 1)), B);
        clear Bindices;
        C=sparse(Cindices(1:2:size(Cindices, 1)-1), Cindices(2:2:size(Cindices, 1)), C);
        clear Cindices;
        invalidEentries=find(Eindices==0);
        Eindices(invalidEentries)=[];
        E(invalidE2entries(2:2:size(invalidEentries, 1))/2)=[];
        E=sparse(Eindices(1:2:size(Eindices, 1)-1), Eindices(2:2:size(Eindices, 1)), E);
        clear Eindices;  
    elseif isfield(opts, 'obsMaskPrecon')
        [B, C, E, B2, C2, E2, g, g2, Bindices, Cindices, Eindices, B2indices, C2indices, E2indices]=calculate_jacobian_croute(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.n_pts, bundleData.cidx, bundleData.pidx, D, opts.obsMaskPrecon, size(opts.obsMaskPrecon, 1), f, f2);
        B2=sparse(B2indices(1:2:size(B2indices, 1)-1), B2indices(2:2:size(B2indices, 1)), B2);
        clear B2indices;
        C2=sparse(C2indices(1:2:size(C2indices, 1)-1), C2indices(2:2:size(C2indices, 1)), C2);
        clear C2indices;
        invalidE2entries=find(E2indices==0);
        E2indices(invalidE2entries)=[];
        E2(invalidE2entries(2:2:size(invalidE2entries, 1))/2)=[];
        E2=sparse(E2indices(1:2:size(E2indices, 1)-1), E2indices(2:2:size(E2indices, 1)), E2);
        clear E2indices;         
        diag_weights=0;
        B=sparse(Bindices(1:2:size(Bindices, 1)-1), Bindices(2:2:size(Bindices, 1)), B);
        clear Bindices;
        C=sparse(Cindices(1:2:size(Cindices, 1)-1), Cindices(2:2:size(Cindices, 1)), C);
        clear Cindices;
        E=sparse(Eindices(1:2:size(Eindices, 1)-1), Eindices(2:2:size(Eindices, 1)), E);
        clear Eindices;   
    else
        [B, C, E, B2, C2, E2, g, g2, Bindices, Cindices, Eindices, ~, ~, ~]=calculate_jacobian_croute(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.n_pts, bundleData.cidx, bundleData.pidx, D, f);
        diag_weights=0;
        B=sparse(Bindices(1:2:size(Bindices, 1)-1), Bindices(2:2:size(Bindices, 1)), B);
        clear Bindices;
        C=sparse(Cindices(1:2:size(Cindices, 1)-1), Cindices(2:2:size(Cindices, 1)), C);
        clear Cindices;
        E=sparse(Eindices(1:2:size(Eindices, 1)-1), Eindices(2:2:size(Eindices, 1)), E);
        clear Eindices;   
    end
%{
    size_p=bundleData.size_p;
    size_c=bundleData.size_c;
    size_o=bundleData.size_o;
    size_x=size_p*bundleData.n_pts+size_c*bundleData.n_cams;
    n_cams=bundleData.n_cams;
    
    g2=zeros(size_x, 1);
    B=zeros(size_c*size_c*bundleData.n_cams, 3);    
    E=zeros(size_c*size_p*bundleData.n_obs, 3);
    C=zeros(size_p*size_p*bundleData.n_pts, 3);    
    for i=1:bundleData.n_obs
        camera=bundleData.cidx(i);
        point=bundleData.pidx(i);
        c=state_x(camera*size_c+1:camera*size_c+9);
        x=state_x(bundleData.n_cams*size_c+point*size_p+1:bundleData.n_cams*size_c+point*size_p+3);                                
        g_buffer=calculate_one_jacobian_backup(c, x, 24);
        
        % What is scaling?
        if size(D, 1)>0
            for j=1:size_o
               for k=1:size_c
                   g_buffer((j-1)*size_c+k)=g_buffer((j-1)*size_c+k)*D(size_c*camera+k);
                end
                for k=1:size_p
                    g_buffer(size_o*size_c+(j-1)*size_p+k)=g_buffer(size_o*size_c+(j-1)*size_p+k)*D(size_c*bundleData.n_cams+size_p*point+k);
                end
            end            
         %   g_buffer(1:size_c)=g_buffer(1:size_c).*D(size_c*camera+1:size_c*camera+size_c);
         %   g_buffer(size_c+1:2*size_c)=g_buffer(size_c+1:2*size_c).*D(size_c*camera+1:size_c*camera+size_c);
         %               
         %   D_indices=size_c*bundleData.n_cams+size_p*point+1:size_c*bundleData.n_cams+size_p*point+size_p;
            
         %   g_indices=size_o*size_c+1:size_o*size_c+size_p;
         %   g_buffer(g_indices)=g_buffer(g_indices).*D(D_indices);
            
         %   g_indices=size_o*size_c+size_p+1:size_o*size_c+size_p*2;
         %   g_buffer(g_indices)=g_buffer(g_indices).*D(D_indices);
        end
      
    %     B(camera*size_c+1:(camera+1)*size_c,camera*size_c+1:(camera+1)*size_c)=g_buffer(1:size_c*size_o)*g_buffer(1:size_c*size_o)'
        for blas_i=1:size_c
            for blas_j=1:size_c             
                for blas_k=1:size_o         
                    B(camera*size_c*size_c+(blas_j-1)*size_c+blas_i, 1:2)=[camera*size_c+blas_i camera*size_c+blas_j];
                    B(camera*size_c*size_c+(blas_j-1)*size_c+blas_i, 3)=B(camera*size_c*size_c+(blas_j-1)*size_c+blas_i, 3)+g_buffer((blas_k-1)*size_c+blas_i)*g_buffer((blas_k-1)*size_c+blas_j);
                end
            end
        end                        
    %    BBlock=g_buffer(1:size_c)*g_buffer(1:size_c)'+g_buffer(size_c+1:size_c*2)*g_buffer(size_c+1:size_c*2)';
    %    [row, col, val]=find(BBlock);
    %    B(camera*size_c*size_c+1:(camera+1)*size_c*size_c, 1:2)=[row+camera*size_c col+camera*size_c];
    %    B(camera*size_c*size_c+1:(camera+1)*size_c*size_c, 3)=B(camera*size_c*size_c+1:(camera+1)*size_c*size_c, 3)+val;
        
        
        for blas_i=1:size_p
            for blas_j=1:size_c
                for blas_k=1:size_o                   
                    E((i-1)*size_p*size_c+(blas_j-1)*size_p+blas_i, 1:2)=[camera*size_c+blas_j point*size_p+blas_i];
                    E((i-1)*size_p*size_c+(blas_j-1)*size_p+blas_i, 3)=E((i-1)*size_p*size_c+(blas_j-1)*size_p+blas_i, 3)+1.0*g_buffer(size_c*size_o+(blas_k-1)*size_p+blas_i)*g_buffer((blas_k-1)*size_c+blas_j);
                end
            end
        end
    %    EBlock=g_buffer(1:size_c)*g_buffer(size_c*2+1:size_c*2+size_p)'+g_buffer(size_c+1:size_c*2)*g_buffer(size_c*2+size_p+1:size_c*2+size_p*2)';
    %    [row, col, val]=find(EBlock);
    %    E((i-1)*size_p*size_c+1:i*size_p*size_c, 1:2)=[row+camera*size_c col+point*size_p];
    %    E((i-1)*size_p*size_c+1:i*size_p*size_c, 3)=E((i-1)*size_p*size_c+1:i*size_p*size_c, 3)+val;
        
   %      C(point*size_p+1:(point+1)*size_p,
   %      point*size_p+1:(point+1)*size_p)=g_buffer(size_c*size_o+1:size_o*(size_c+size_p))*g_buffer(size_c*size_o+1:size_o*(size_c+size_p))'
        for blas_i=1:size_p
            for blas_j=1:size_p
                for blas_k=1:size_o
                    C(point*size_p*size_p+(blas_j-1)*size_p+blas_i, 1:2)=[point*size_p+blas_i point*size_p+blas_j];
                    C(point*size_p*size_p+(blas_j-1)*size_p+blas_i, 3)=C(point*size_p*size_p+(blas_j-1)*size_p+blas_i, 3)+g_buffer(size_c*size_o+(blas_k-1)*size_p+blas_i)*g_buffer(size_c*size_o+(blas_k-1)*size_p+blas_j);
                end
            end
        end
    %    CBlock=g_buffer(size_c*2+1:size_c*2+size_p)*g_buffer(size_c*2+1:size_c*2+size_p)'+g_buffer(size_c*2+size_p+1:size_c*2+size_p*2)*g_buffer(size_c*2+size_p+1:size_c*2+size_p*2)';
    %    [row, col, val]=find(CBlock);
    %    C(point*size_p*size_p+1:(point+1)*size_p*size_p, 1:2)=[row+point*size_p col+point*size_p];
    %    C(point*size_p*size_p+1:(point+1)*size_p*size_p, 3)=C(point*size_p*size_p+1:(point+1)*size_p*size_p, 3)+val;
        
        for blas_i=1:size_c
            for blas_j=1:size_o
                g2(camera*size_c+blas_i)=g2(camera*size_c+blas_i)-g_buffer((blas_j-1)*size_c+blas_i)*f((i-1)*size_o+blas_j);
            end
        
        %   above is equivalent to g(camera*size_c+blas_i)=g(camera*size_c+blas_i)-g_buffer(blas_i)*f(i*2-1)-g_buffer(size_c+blas_i)*f(i*2);
        end
    %    g2(camera*size_c+1:camera*size_c+size_c)=g2(camera*size_c+1:camera*size_c+size_c)-f(i*2-1)*g_buffer(1:size_c)-f(i*2)*g_buffer(size_c+1:size_c*2);
        
        for blas_i=1:size_p
            for blas_j=1:size_o
                g2(bundleData.n_cams*size_c+point*size_p+blas_i)=g2(bundleData.n_cams*size_c+point*size_p+blas_i)-g_buffer(size_c*size_o+(blas_j-1)*size_p+blas_i)*f((i-1)*size_o+blas_j);
            end
       %     above is equivalent to
       %     g(bundleData.n_cams*size_c+point*size_p+blas_i)=g(bundleData.n_cams*size_c+point*size_p+blas_i)-g_buffer(size_c*size_o+blas_i)*f((i*2-1)-g_buffer(size_c*size_o+size_p+blas_i)*f(i*2);
        end
    %    g2(n_cams*size_c+point*size_p+1:n_cams*size_c+point*size_p+size_p)=g2(n_cams*size_c+point*size_p+1:n_cams*size_c+point*size_p+size_p)-f(i*2-1)*g_buffer(size_c*2+1:size_c*2+size_p)-f(i*2)*g_buffer(size_c*2+size_p+1:size_c*2+size_p*2);
    end
    B=sparse(B(:, 1), B(:, 2), B(:, 3));
    C=sparse(C(:, 1), C(:, 2), C(:, 3));
    E=sparse(E(:, 1), E(:, 2), E(:, 3));
%}
end