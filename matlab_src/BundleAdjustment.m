classdef BundleAdjustment
    properties
        bundleData
        
        B
        C
        Ci
        E
        D
        
        f
        g
        
        S
        b
        newS
        
        newCamMat
        diag_weights
        
        state_x
        
        % LM iteration related parameters
        LM_max_iterCount
        tau
        tolFun
        tolX
        min_rho
        mu
        mu_multiplier
        
        old_pcg_solution
        pcgtol
        maxPCGIter
        
        opts
    end
    
    methods
        function obj=BundleAdjustment(opts)
            if isfield(opts, 'bundleData')
                bundle_data=opts.bundleData;
            else
                bundle_data=read_bundle_data(opts);    
            end
            obj.state_x=bundle_data.x;
            clear bundle_data.x;
            obj.LM_max_iterCount=10;                        
            obj.bundleData=bundle_data;
            obj.tau=1e-4;
            obj.tolFun=1e-10;
            obj.tolX=1e-10;
            obj.min_rho=1e-3;
            obj.mu_multiplier=2;
            
            obj.pcgtol=opts.pcgtol;
            obj.maxPCGIter=500;
            obj.old_pcg_solution=zeros(bundle_data.n_cams*9, 1);
            obj.opts=opts;
        end
        
        function initLM(obj)   
            size_c=obj.bundleData.size_c;
            size_p=obj.bundleData.size_p;
            size_o=obj.bundleData.size_o;    
            size_x=size_p*obj.bundleData.n_pts + size_c*obj.bundleData.n_cams;

            obj.D=ones(size_x, 1);
            pt_cam=cell(obj.bundleData.n_pts, 1);
            for i=1:obj.bundleData.n_obs
                camera=obj.bundleData.cidx(i);
                point=obj.bundleData.pidx(i);
                pt_cam{point+1}=[pt_cam{point+1}; camera+1];

                if obj.opts.use_scaling
                    c=obj.state_x(camera*size_c+1:camera*size_c+9);
                    x=obj.state_x(obj.bundleData.n_cams*size_c+point*size_p+1:obj.bundleData.n_cams*size_c+point*size_p+3);
                    g_buffer=calculate_one_jacobian(c, x);                
                    for j=1:size_o
                        for k=1:size_c
                            obj.D(size_c*camera+k)=obj.D(size_c*camera+k)+g_buffer((j-1)*size_c+k)^2;
                        end

                        for k=1:size_p
                            obj.D(size_c*obj.bundleData.n_cams+size_p*point+k)=obj.D(size_c*obj.bundleData.n_cams+size_p*point+k)+g_buffer(size_c*size_o+(j-1)*size_p+k)^2;
                        end
                    end
                end
            end

            obj.newCamMat=build_camera_matrix(pt_cam, obj.bundleData.n_cams, obj.opts);
            if obj.opts.use_scaling
                for i=1:size_x
                    obj.D(i)=1/(1e-18+sqrt(obj.D(i)));
                end
            end
        end

        function solution_status=LM(obj)                       
            solution_status=0;    
            
            if isfield(obj.opts, 'newCamMat') && isfield(obj.opts, 'diagonalScale');
                obj.D=obj.opts.diagonalScale;
                obj.newCamMat=obj.opts.newCamMat;
                clear opts.newCamMat;
                clear opts.diagonalScale;
            else
                initLM(obj);
            end
                    
            obj.f=calculate_residual(obj.state_x, obj.bundleData.n_cams, obj.bundleData.n_obs, obj.bundleData.cidx, obj.bundleData.pidx, obj.bundleData.obs);
            [obj.B, obj.C, obj.E, obj.g, obj.diag_weights]=calculate_jacobian(obj.state_x, obj.bundleData, obj.f, obj.D, obj.opts);

            objval=obj.f'*obj.f;
            gnorm=norm(obj.g, 'inf');
            xnorm=norm(obj.state_x, 2);

            LM_converged=(gnorm<obj.tolFun);    
            if LM_converged
                solution_status=3;
            end

            obj.mu=obj.tau;
            k=0;     

            if obj.opts.dumpSolutionStatus==1
                fprintf(obj.opts.dumpFile, '% 4d: f:% 8e g:% 8e mu:% 8e \n',k,objval, gnorm, obj.mu);
            else
                fprintf('% 4d: f:% 8e g:% 8e mu:% 8e \n',k,objval, gnorm, obj.mu);
            end
            
            while (~LM_converged)&&(k<obj.LM_max_iterCount)
                k=k+1;
                LM_step_trials=0;
                [LM_step_converged, h, schur_model_reduction]=obj.calculate_LM_step();
                while (~LM_step_converged)&&(LM_step_trials<=10)
                    obj.mu=obj.mu*obj.mu_multiplier;
                    obj.mu_multiplier=obj.mu_multiplier*2;
                    LM_step_trials=LM_step_trials+1;
                    [LM_step_converged, h, schur_model_reduction]=obj.calculate_LM_step();
                end

                if(LM_step_trials==11)
                    solution_status=4; % calculate_LM_step failed to find a good step size
                    return;
                end

                hd=h.*obj.D;
                hnorm=norm(hd, 2);
                state_xn=obj.state_x+hd;        

                if hnorm<=obj.tolX*(xnorm+obj.tolX)
                    LM_converged=true;
                    solution_status=2; % Solution found. Terminated by tolX
                end
             %   fn2=calculate_residual_backup(state_xn, bundleData);
                fn=calculate_residual(state_xn, obj.bundleData.n_cams, obj.bundleData.n_obs, obj.bundleData.cidx, obj.bundleData.pidx, obj.bundleData.obs);
                objval_n=norm(fn, 2)^2;
                rho=(objval-objval_n)/schur_model_reduction;
                if rho>obj.min_rho
                    obj.state_x=state_xn;
                    objval=objval_n;
                    obj.f=fn;                   
                    xnorm=norm(obj.state_x, 2);
                    [obj.B, obj.C, obj.E, obj.g, obj.diag_weights]=calculate_jacobian(obj.state_x, obj.bundleData, obj.f, obj.D);
                    gnorm=norm(obj.g, 'inf');
                    if gnorm<=obj.tolFun
                        LM_converged=true;
                        solution_status=3; % Solution found. Terminated by tolFun
                    end
                    obj.mu=obj.mu*max(1.0/3, 1-power(2*rho-1, 3));
                    obj.mu_multiplier=2;
                else
                    obj.mu=obj.mu*obj.mu_multiplier;
                    obj.mu_multiplier=obj.mu_multiplier*2;
                end

                if obj.opts.dumpSolutionStatus==1
                    fprintf(obj.opts.dumpFile, '% 4d: f:% 8e g:% 8e h:% 8e rho:% 8e mu:% 8e \n',k,objval, gnorm,hnorm, rho, obj.mu);
                else
                    fprintf('% 4d: f:% 8e g:% 8e h:% 8e rho:% 8e mu:% 8e \n',k,objval, gnorm,hnorm, rho, obj.mu);
                end

                if k>=obj.LM_max_iterCount
                    solution_status=1; % Max iterations for LM reached. Iteration terminated
                end
            end
        end
        
        function [LM_step_converged, h, schur_model_reduction, delta_y]=calculate_LM_step(obj)               
            n_cams=obj.bundleData.n_cams;
            size_c=obj.bundleData.size_c;

            % Scale the diagonal of the Jacobian, i.e B and C. 
            if ~obj.opts.largeScale
                min_diag=1e-6;
                max_diag=1e+32;
                obj.diag_weights=obj.mu*max(min(diag(obj.B(1:9, 1:9)), max_diag), min_diag);
                obj.B=obj.B+diag(obj.mu*max(min(diag(obj.B), max_diag), min_diag));
                obj.C=obj.C+diag(obj.mu*max(min(diag(obj.C), max_diag), min_diag));
            end    

            % Ci is the block-inverse of C   
            if ~obj.opts.largeScale
                [Crow, Ccol, Cval]=find(obj.C);
                Cval=C_diagonal_inverse(Cval);
                obj.Ci=sparse(Crow, Ccol, Cval);
                clear Crow;clear Ccol;clear Cval;
                ECi=obj.E*obj.Ci;
                obj.S=obj.B-ECi*obj.E';        
                obj.b=obj.g(1:n_cams*size_c)-ECi*obj.g(n_cams*size_c+1:size(obj.g, 1));

                [newCamMat_row, newCamMat_col]=find(obj.newCamMat);        
                obj.newS=spalloc(size(obj.S, 1), size(obj.S, 2), size(newCamMat_row, 1)*size_c*size_c);
                for i=1:size(newCamMat_row, 1)
                    rowBlock=(newCamMat_row(i)-1)*size_c+1:newCamMat_row(i)*size_c;
                    colBlock=(newCamMat_col(i)-1)*size_c+1:newCamMat_col(i)*size_c;
                    obj.newS(rowBlock, colBlock)=obj.S(rowBlock, colBlock);
                end
                clear newCamMat_row;clear newCamMat_col;        
            else
                obj.Ci=C_diagonal_inverse(obj.C);
                obj.S=S_handle_large_scale();
                obj.b=build_b_large_scale(obj.E, obj.Ci, obj.g, obj.bundleData.n_cams, obj.bundleData.n_obs, obj.bundleData.n_pts, obj.bundleData.cidx, obj.bundleData.pidx);
                obj.newS=build_preconditioner_large_scale();
            end

            if obj.opts.reducedCameraMatrixUsage==0 % Do not use newS
                if obj.opts.solverMethod==1 % PCG
                    if obj.opts.PCGUsePreviousSolution
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.S, obj.b, obj.pcgtol, obj.maxPCGIter, [], [], obj.old_pcg_solution);
                    else    
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.S, obj.b, obj.pcgtol, obj.maxPCGIter);
                    end
                elseif obj.opts.solverMethod==2 % Direct solve
                    delta_y=obj.S\obj.b;
                end
            elseif obj.opts.reducedCameraMatrixUsage==1 % newS as preconditioner
                if obj.opts.solverMethod==1 % PCG
                    if obj.opts.PCGUsePreviousSolution
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.S, obj.b, obj.pcgtol, obj.maxPCGIter, obj.newS, [], obj.old_pcg_solution);
                    else                
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.S, obj.b, obj.pcgtol, obj.maxPCGIter, obj.newS);
                    end
                elseif obj.opts.solverMethod==2 % Direct solve
                    delta_y=obj.S\obj.b;
                end
            elseif obj.opts.reducedCameraMatrixUsage==2 % newS as solver
                if obj.opts.solverMethod==1 % PCG
                    if obj.opts.PCGUsePreviousSolution
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.newS, obj.b, obj.pcgtol, obj.maxPCGIter, [], [], obj.old_pcg_solution);
                    else
                        [delta_y, pcg_flag, pcg_res, pcg_iter, normrs]=my_pcg(obj.newS, obj.b, obj.pcgtol, obj.maxPCGIter);
                    end
                elseif obj.opts.solverMethod==2 % Direct solve
                    delta_y=obj.newS\obj.b;
                end
            end
            obj.old_pcg_solution=delta_y;
            
            fprintf(obj.opts.normrsDump, '%d\n', size(normrs, 1));
            fprintf(obj.opts.normrsDump, '%d\n', normrs);

            if obj.opts.reducePCGTol        
                obj.pcgtol=obj.pcgtol/2;
            end

            if obj.opts.solverMethod==1
                if obj.opts.dumpSolutionStatus==1 
                    fprintf(obj.opts.dumpFile, 'pcg flag: %d pcg residual: %d pcg iteration: %d\n', pcg_flag, pcg_res, pcg_iter);
                else
                    fprintf('pcg flag: %d pcg residual: %d pcg iteration: %d\n', pcg_flag, pcg_res, pcg_iter);
                end
            end

            delta_z=obj.Ci*obj.g(n_cams*size_c+1:size(obj.g, 1))-ECi'*delta_y;
            clear ECi;
            h=[delta_y; delta_z];

            schur_model_reduction=obj.g(n_cams*size_c+1:size(obj.g, 1))'*obj.Ci*obj.g(n_cams*size_c+1:size(obj.g, 1));
            clear Ci;
            if obj.opts.reducedCameraMatrixUsage==0 || obj.opts.reducedCameraMatrixUsage==1
                schur_model_reduction=schur_model_reduction+(obj.b*2-obj.S*delta_y)'*delta_y;
            elseif obj.opts.reducedCameraMatrixUsage==2
                schur_model_reduction=schur_model_reduction+(obj.b*2-obj.newS*delta_y)'*delta_y;
            end
            schur_model_reduction=schur_model_reduction+(h(1:9).*h(1:9))'*obj.diag_weights(1:9);
            LM_step_converged=true;
        end
        
        % Computes and returns y=Ax, where A=B-E'CiE
        function y=eval_Axb(obj, x)
            y=eval_Axb_croute();
        end
    end   
end