function bundlerBatchjob
    opts=struct('data_filename', 'problem-3068-310854-pre.txt', 'reducedCameraMatrixType', 1, 'reducedCameraMatrixUsage', 1, 'solverMethod', 1, 'dumpSolutionStatus', 1, 'spanner_t', 13, 'use_scaling', true);opts.dumpFile=fopen('4585_1324582-CamMat1-CamUsg1-Spt13.txt', 'w');
    runBundler(opts); 
    opts=struct('data_filename', 'problem-3068-310854-pre.txt', 'reducedCameraMatrixType', 2, 'reducedCameraMatrixUsage', 1, 'solverMethod', 1, 'dumpSolutionStatus', 1, 'spanner_t', 13, 'use_scaling', true);opts.dumpFile=fopen('4585_1324582-CamMat2-CamUsg1-Spt13.txt', 'w');
    runBundler(opts);
    opts=struct('data_filename', 'problem-3068-310854-pre.txt', 'reducedCameraMatrixType', 1, 'reducedCameraMatrixUsage', 2, 'solverMethod', 1, 'dumpSolutionStatus', 1, 'spanner_t', 13, 'use_scaling', true);opts.dumpFile=fopen('4585_1324582-CamMat1-CamUsg2-Spt13.txt', 'w');
    runBundler(opts);
    opts=struct('data_filename', 'problem-3068-310854-pre.txt', 'reducedCameraMatrixType', 3, 'reducedCameraMatrixUsage', 1, 'solverMethod', 1, 'dumpSolutionStatus', 1, 'spanner_t', 16, 'use_scaling', true);opts.dumpFile=fopen('4585_1324582-CamMat1-CamUsg2-Spt16.txt', 'w');
    runBundler(opts);
end