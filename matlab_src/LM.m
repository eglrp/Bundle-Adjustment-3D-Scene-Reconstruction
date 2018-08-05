function [state_x, solution_status]=LM(bundleData, camMat, opts)
    state_x=bundleData.x;
    clear bundleData.x;
    tau=1e-4;
    tolFun=1e-10;
    tolX=1e-10;
    min_rho=1e-3;    
    mu_multiplier=2;
    solution_status=0;
    newCamMat=zeros(1);
    if isfield(opts, 'LM_maxIter_count')
        LM_maxIter_count=opts.LM_maxIter_count;
    else
        LM_maxIter_count=10;
    end
    
    if isfield(opts, 'preconMatVarFile') && isfield(opts, 'skeletalMatVarFile')
        var=load(opts.preconMatVarFile);
        var2=load(opts.skeletalMatVarFile);
        opts.obsMask=var2.obs_mask;
        opts.obsMaskPrecon=var.obs_mask;
        D=var2.D;
        newCamMat=var.newCamMat;
    elseif isfield(opts, 'preconMatVarFile')
         var=load(opts.preconMatVarFile);         
         opts.obsMaskPrecon=var.obs_mask;
         newCamMat=var.newCamMat;
         D=var.D;
    else
        if opts.reducedCameraMatrixType>0
            newCamMat=build_camera_matrix(camMat, opts);
        end
        [D, bundleData]=initLM(bundleData, state_x);           
    end
    
    
    if isfield(opts, 'obsMaskPrecon') && isfield(opts, 'obsMask')
        [f, f2]=calculate_residual(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs, opts.obsMaskPrecon, opts.obsMask, size(opts.obsMaskPrecon, 1), size(opts.obsMask, 1));
        [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
    elseif isfield(opts, 'obsMaskPrecon')
        [f, f2]=calculate_residual(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs, opts.obsMaskPrecon, size(opts.obsMaskPrecon, 1));
        [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
    else
        [f, f2]=calculate_residual(state_x, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs);
        [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
    end
    
    objval=f'*f;
    gnorm=norm(g, 'inf');
    xnorm=norm(state_x, 2);
    
    LM_converged=(gnorm<tolFun);    
    if LM_converged
        solution_status=3;
    end
        
    mu=tau;
    k=0;     
    
    if isfield(opts, 'dumpSolutionToFile')
        fprintf(opts.dumpSolutionToFile, '%d: f:%e g:%e mu:%e \n',k,objval, gnorm, mu);
    else
        fprintf('%d: f:%e g:%e mu:%e \n',k,objval, gnorm, mu);
    end
    
    delta_y=zeros(bundleData.n_cams*bundleData.size_c, 1);   
    while (~LM_converged)&&(k<LM_maxIter_count)
        k=k+1;
        LM_step_trials=0;
        if isfield(opts, 'obsMaskPrecon')
            [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(bundleData, B, C, E, g, mu, newCamMat, delta_y, diag_weights, opts, B2, C2, E2, g2);
        else
            [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(bundleData, B, C, E, g, mu, newCamMat, delta_y, diag_weights, opts);
        end
        
        while (~LM_step_converged)&&(LM_step_trials<=10)
            mu=mu*mu_multiplier;
            mu_multiplier=mu_multiplier*2;
            LM_step_trials=LM_step_trials+1;
            if isfield(opts, 'obsMaskPrecon')
                [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(bundleData, B, C, E, g, mu, newCamMat, delta_y, diag_weights, opts, B2, C2, E2, g2);
            else
                [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(bundleData, B, C, E, g, mu, newCamMat, delta_y, diag_weights, opts);
            end
        end
        
        if(LM_step_trials==11)
            solution_status=4; % calculate_LM_step failed to find a good step size
            return;
        end
       
        hd=h.*D;
        hnorm=norm(hd, 2);
        state_xn=state_x+hd;        
        
        if hnorm<=tolX*(xnorm+tolX) && isfield(opts, 'earlyTermination')
            LM_converged=true;
            solution_status=2; % Solution found. Terminated by tolX
        end
        
        switchObjective=false;
        if  k==LM_maxIter_count/2 && isfield(opts, 'obsMask')
            opts=rmfield(opts, 'obsMask');
            switchObjective=true;
        end
    
        if isfield(opts, 'obsMaskPrecon') && isfield(opts, 'obsMask')
            [fn, f2n]=calculate_residual(state_xn, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs, opts.obsMaskPrecon, opts.obsMask, size(opts.obsMaskPrecon, 1), size(opts.obsMask, 1));
        elseif isfield(opts, 'obsMaskPrecon')
            [fn, f2n]=calculate_residual(state_xn, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs, opts.obsMaskPrecon, size(opts.obsMaskPrecon, 1));
        else
            [fn, f2n]=calculate_residual(state_xn, bundleData.n_cams, bundleData.n_obs, bundleData.cidx, bundleData.pidx, bundleData.obs);
        end
        objval_n=norm(fn, 2)^2;
        rho=(objval-objval_n)/schur_model_reduction;
        if rho>min_rho || switchObjective
            state_x=state_xn;
            objval=objval_n;
            f=fn;
            f2=f2n;
            xnorm=norm(state_x, 2);            
            gnorm=norm(g, 'inf');
            
            if isfield(opts, 'obsMaskPrecon') && isfield(opts, 'obsMask')
                [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
            elseif isfield(opts, 'obsMaskPrecon')
                [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
            else
                [B, C, E, g, B2, E2, C2, g2, diag_weights]=calculate_jacobian(state_x, bundleData, f, f2, D, opts);
            end
            if gnorm<=tolFun && isfield(opts, 'earlyTermination')
                LM_converged=true;
                solution_status=3; % Solution found. Terminated by tolFun
            end
            if switchObjective
                % Don't change mu
            else
                mu=mu*max(1.0/3, 1-power(2*rho-1, 3));
            end
            mu_multiplier=2;
        else
            mu=mu*mu_multiplier;
            mu_multiplier=mu_multiplier*2;
        end
        
        if isfield(opts, 'camDump')
            fprintf(opts.camDump, '%e\n', state_x(1:bundleData.n_cams*9));
        end
        
        if isfield(opts, 'dumpSolutionToFile')
            fprintf(opts.dumpSolutionToFile, '%d: f:%e g:%e h:%e rho:%e mu:%e \n',k,objval, gnorm,hnorm, rho, mu);
        else
            fprintf('%d: f:%e g:%e h:%e rho:%e mu:%e \n',k,objval, gnorm,hnorm, rho, mu);
        end
        
        if k>=LM_maxIter_count
            solution_status=1; % Max iterations for LM reached. Iteration terminated
        end
    end
end

% Jacobian*h=-g
function [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(bundleData, B, C, E, g, mu, newCamMat, delta_y, diag_weights, opts, B2, C2, E2, g2)               
    n_cams=bundleData.n_cams;
    size_c=bundleData.size_c;

    % Scale the diagonal of the Jacobian, i.e B and C. 
    min_diag=1e-6;
    max_diag=1e+32;
    diag_weights=mu*max(min(diag(B(1:9, 1:9)), max_diag), min_diag);
    B=B+diag(mu*max(min(diag(B), max_diag), min_diag));
    C=C+diag(mu*max(min(diag(C), max_diag), min_diag));

    if nargin==15          
        B2=B2+diag(mu*max(min(diag(B2), max_diag), min_diag));        
        C2=C2+diag(mu*max(min(diag(C2), max_diag), min_diag));
    end

    % Ci is the block-inverse of C   
    [Crow, Ccol, Cval]=find(C);
    Cval=C_diagonal_inverse(Cval);
    Ci=sparse(Crow, Ccol, Cval);
    clear C;clear Crow;clear Ccol;clear Cval;
    ECi=E*Ci;
    S=B-ECi*E';
    clear B;        
    b=g(1:n_cams*size_c)-ECi*g(n_cams*size_c+1:size(g, 1));
    if nargin==15
        [C2row, C2col, C2val]=find(C2);
        C2val=C_diagonal_inverse(C2val);
        Ci2=sparse(C2row, C2col, C2val);
        clear C2;clear C2row;clear C2col;clear C2val;
        S2=B2-E2*Ci2*E2';
        clear B;clear Ci2;                

        if opts.reducedCameraMatrixType>0
            [newCamMat_row, newCamMat_col]=find(newCamMat);        
            newS=spalloc(size(S2, 1), size(S2, 2), size(newCamMat_row, 1)*size_c*size_c);
            for i=1:size(newCamMat_row, 1)
                rowBlock=(newCamMat_row(i)-1)*size_c+1:newCamMat_row(i)*size_c;
                colBlock=(newCamMat_col(i)-1)*size_c+1:newCamMat_col(i)*size_c;
                newS(rowBlock, colBlock)=S2(rowBlock, colBlock); %#ok<SPRIX>
            end
            clear newCamMat_row;clear newCamMat_col; 
        end
    else
        if opts.reducedCameraMatrixType>0
            [newCamMat_row, newCamMat_col]=find(newCamMat);        
            newS=spalloc(size(S, 1), size(S, 2), size(newCamMat_row, 1)*size_c*size_c);
            for i=1:size(newCamMat_row, 1)
                rowBlock=(newCamMat_row(i)-1)*size_c+1:newCamMat_row(i)*size_c;
                colBlock=(newCamMat_col(i)-1)*size_c+1:newCamMat_col(i)*size_c;
                newS(rowBlock, colBlock)=S(rowBlock, colBlock); %#ok<SPRIX>
            end
            clear newCamMat_row;clear newCamMat_col;  
        end
    end 
    
    if opts.reducedCameraMatrixUsage==0 % Do not use newS, solve Sx=b directly
        delta_y=S\b;    
    elseif opts.reducedCameraMatrixUsage==1 % Do not use newS, solve Sx=b via PCG
        [delta_y, pcg_flag, pcg_res, pcg_iter]=my_pcg(S, b, opts.pcgtol, 500);
    elseif opts.reducedCameraMatrixUsage==2 % newS as preconditioner
        [delta_y, pcg_flag, pcg_res, pcg_iter]=my_pcg(S, b, opts.pcgtol, 500, newS);        
    elseif opts.reducedCameraMatrixUsage==3 % solve newSx=b via PCG
        [delta_y, pcg_flag, pcg_res, pcg_iter]=my_pcg(newS, b, opts.pcgtol, 500);  
    elseif opts.reducedCameraMatrixUsage==4 % solve newSx=b directly
        delta_y=newS\b;
    end
    
    if opts.reducedCameraMatrixUsage>=1 && opts.reducedCameraMatrixUsage<=3
        if isfield(opts, 'dumpSolutionToFile')
            fprintf(opts.dumpSolutionToFile, 'pcg flag: %d pcg residual: %e pcg iteration: %d\n', pcg_flag, pcg_res, pcg_iter);
        else
            fprintf('pcg flag: %d pcg residual: %e pcg iteration: %d\n', pcg_flag, pcg_res, pcg_iter);
        end
    end
    
    delta_z=Ci*g(n_cams*size_c+1:size(g, 1))-ECi'*delta_y;
    clear E;clear ECi;
    h=[delta_y; delta_z];                                   
    
    schur_model_reduction=g(n_cams*size_c+1:size(g, 1))'*Ci*g(n_cams*size_c+1:size(g, 1));
    clear Ci;
    if opts.reducedCameraMatrixUsage>=0 && opts.reducedCameraMatrixUsage<=2
        schur_model_reduction=schur_model_reduction+(b*2-S*delta_y)'*delta_y;
    elseif opts.reducedCameraMatrixUsage>=3 && opts.reducedCameraMatrixUsage<=4
        schur_model_reduction=schur_model_reduction+(b*2-newS*delta_y)'*delta_y;
    end
    schur_model_reduction=schur_model_reduction+(h(1:9).*h(1:9))'*diag_weights(1:9);
    LM_step_converged=true;
end
    
function [D, bundleData]=initLM(bundleData, state_x)    
    size_c=bundleData.size_c;
    size_p=bundleData.size_p;
    size_o=bundleData.size_o;    
    size_x=size_p*bundleData.n_pts + size_c*bundleData.n_cams;
    pidx=bundleData.pidx;
    cidx=bundleData.cidx;
    
    D=ones(size_x, 1);
    for i=1:bundleData.n_obs
        camera=cidx(i);
        point=pidx(i);

        c=state_x(camera*size_c+1:camera*size_c+9);
        x=state_x(bundleData.n_cams*size_c+point*size_p+1:bundleData.n_cams*size_c+point*size_p+3);
        g_buffer=calculate_one_jacobian(c, x);                
        for j=1:size_o
            for k=1:size_c
                D(size_c*camera+k)=D(size_c*camera+k)+g_buffer((j-1)*size_c+k)^2;
            end

            for k=1:size_p
                D(size_c*bundleData.n_cams+size_p*point+k)=D(size_c*bundleData.n_cams+size_p*point+k)+g_buffer(size_c*size_o+(j-1)*size_p+k)^2;
            end
        end    
    end
  
    for i=1:size_x
        D(i)=1/(1e-18+sqrt(D(i)));
    end
end

