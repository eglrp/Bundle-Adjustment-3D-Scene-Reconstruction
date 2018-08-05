function build_skeletal_on_obs_graph(data_filename, output_filename, t)
    bundleData=struct('size_c', 9, 'size_p', 3, 'size_o', 2);
    bundleData.n_cams=0;
    bundleData.n_pts=0;
    bundleData.n_obs=0;
    fid=fopen(data_filename, 'r');    
    bundleData.n_cams=fscanf(fid, '%d', 1);
    bundleData.n_pts=fscanf(fid, '%d', 1);
    bundleData.n_obs=fscanf(fid, '%d', 1);
    
    S=spalloc(bundleData.n_cams+bundleData.n_pts, bundleData.n_cams+bundleData.n_pts, bundleData.n_obs*2);
    for i=1:bundleData.n_obs
        cid=fscanf(fid, '%d', 1);
        pid=fscanf(fid, '%d', 1);
        obs_in1=fscanf(fid, '%f', 1);
        obs_in2=fscanf(fid, '%f', 1);
        S(cid+1, bundleData.n_cams+pid+1)=1;
    end
    S=S+S';
    opts=struct('reducedCameraMatrixType', 3, 'spanner_t', t);
    G=build_camera_matrix(S, opts);
    save(output_filename, 'G');
    fclose(fid);
end