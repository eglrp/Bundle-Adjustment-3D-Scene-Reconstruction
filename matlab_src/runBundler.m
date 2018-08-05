%{ 
 opts fields:
 
 ======required fields======

 data_filename: filename for reading in the data

 reducedCameraMatrixType: 
 0=don't build reduced camera matrix
 1=diagonal (block jacobi)
 2=max-leaf spanning tree
 3=skeletal graph

 reducedCameraMatrixUsage: 
 0=do not use reduced camera matrix, solve system with \
 1=do not use reduced camera matrix, solve system with PCG
 2=use reduced camera matrix as preconditioner and solve with PCG
 3=use reduced camera matrix to solve the system with PCG
 4=use reduced camera matrix to solve the system with \

 spanner_t: the stretch factor used in building a skeletal graph (only
 required when reducedCameraMatrixType==3)

 pcgtol: tolerance for PCG (only required when
 1<=reducedCameraMatrixType<=3)

 ======optional fields======

 LM_max_iterCount: Maximum LM iteration count

 dumpSolutionToFile: if exists, dump solution and all output to file named
 by this field
 
 camDump: if exists, dumps the camera parameters to this file at every LM
 iteration

 earlyTermination: if exists, terminate LM if it satisfy certain
 inequalities; otherwise, run LM till LM_max_iterCount is reached

 dropObservationsFromFile: if exists, drop observations based on a
 skeletalized bipartite observation graph

 addObsBack: adds back observations and rerun LM. addObsBack stores the LM
 iterations after observations have been added back

 preconMatVarFile, skeletalMatVarFile: currently unused, but the code is
 still there in case they're used in the future. DO NOT add these two
 fields to opts, otherwises things may mess up.
%}
function x=runBundler(opts)        
    [bundleData, camMat]=read_bundle_data(opts);     
    [x, ~]=LM(bundleData, camMat, opts);     
    if isfield(opts, 'addObsBack')
        opts.LM_max_iterCount=opts.addObsBack;
        opts=rmfield(opts, 'addObsBack');
        if isfield(opts, 'dropObservationsFromFile')
            opts=rmfield(opts, 'dropObservationsFromFile');
        end
        opts.reducedCameraMatrixType=2;
        opts.reducedCameraMatrixUsage=2;
        opts.pcgtol=0.1;
        clear bundleData camMat;
        [bundleData, camMat]=read_bundle_data(opts);
        bundleData.x=x;
        clear x;
        [x, ~]=LM(bundleData, camMat, opts);
    end
    if isfield(opts, 'dumpSolutionToFile')
        fprintf(opts.dumpSolutionToFile, '\n\n\n\n\n');
        fprintf(opts.dumpSolutionToFile, '%d\n', x(:));
    end
 end

