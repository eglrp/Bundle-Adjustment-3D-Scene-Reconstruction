% Notes:
% 1: Code doesn't support prior_mask, prior_value, prior_scale and fix_mask yet
% 2: use_scaling is always enabled
% 3: Currently, code would only work for size_c=9, size_p=3, size_o=2
function [bundleData, camMat]=read_bundle_data(opts)
    data_filename=opts.data_filename;
    bundleData=struct('size_c', 9, 'size_p', 3, 'size_o', 2);
    bundleData.n_cams=0;
    bundleData.n_pts=0;
    bundleData.n_obs=0;
    fid=fopen(data_filename, 'r');    
    bundleData.n_cams=fscanf(fid, '%d', 1);
    bundleData.n_pts=fscanf(fid, '%d', 1);
    bundleData.n_obs=fscanf(fid, '%d', 1);
    bundleData.pidx=zeros(bundleData.n_obs, 1);
    bundleData.cidx=-ones(bundleData.n_obs, 1);
    bundleData.obs=zeros(bundleData.n_obs*2, 1);
    size_x=bundleData.n_cams*9+bundleData.n_pts*3;
    bundleData.x=zeros(size_x, 1);
    
    if isfield(opts, 'dropObservationsFromFile')        
        fid2=fopen(opts.dropObservationsFromFile, 'r');        
        cidx=zeros(bundleData.n_obs, 1);pidx=zeros(bundleData.n_obs, 1);
        i=1;      
        while ~feof(fid2)
            cid=fscanf(fid2, '%d', 1);
            pid=fscanf(fid2, '%d', 1);
            if size(cid, 1)==0 || size(pid, 1)==0
                break;
            end
            cidx(i)=cid+1;pidx(i)=pid+1;  
            i=i+1;
        end
        cidx(i:bundleData.n_obs)=[];
        pidx(i:bundleData.n_obs)=[];
        fclose(fid2);
        obsGraph=sparse(cidx, pidx, ones(size(cidx, 1), 1));
        clear cidx pidx;
    end
    
    pointTracks=cell(bundleData.n_pts, 1);
    obs_count=1;
    for i=1:bundleData.n_obs
        cid=fscanf(fid, '%d', 1);
        pid=fscanf(fid, '%d', 1);
        obs_in1=fscanf(fid, '%f', 1);
        obs_in2=fscanf(fid, '%f', 1);
    %    obs_in1=fscanf(fid, '%f', 1)+11*(0.5-rand(1));
    %    obs_in2=fscanf(fid, '%f', 1)+11*(0.5-rand(1));
        if isfield(opts, 'dropObservationsFromFile')
            if obsGraph(cid+1, pid+1)>0
                bundleData.pidx(obs_count)=pid;
                bundleData.cidx(obs_count)=cid;
                bundleData.obs(obs_count*2-1:obs_count*2)=[obs_in1; obs_in2]; 
                pointTracks{pid+1}=[pointTracks{pid+1}; cid+1 i];
                obs_count=obs_count+1;
            end
        else
            bundleData.pidx(i)=pid;
            bundleData.cidx(i)=cid;
            bundleData.obs(i*2-1:i*2)=[obs_in1; obs_in2]; 
            pointTracks{pid+1}=[pointTracks{pid+1}; cid+1 i];
        end
    end    
    
    if isfield(opts, 'state_x')
        bundleData.x=opts.state_x;
    else
        for i=1:size_x        
            x_in=fscanf(fid, '%f', 1);        
            bundleData.x(i)=x_in;  
        end 
    end
    
    if isfield(opts, 'dropObservationsFromFile')        
        bundleData.cidx(obs_count:bundleData.n_obs)=[];
        bundleData.pidx(obs_count:bundleData.n_obs)=[];
        bundleData.obs([obs_count*2-1:bundleData.n_obs*2])=[];
        bundleData.n_obs=obs_count-1;
    end 
    [bundleData.pidx sort_idx]=sort(bundleData.pidx);
    bundleData.cidx=bundleData.cidx(sort_idx);
    bundleData.obs(1:2:bundleData.n_obs*2-1)=bundleData.obs(sort_idx*2-1);
    bundleData.obs(2:2:bundleData.n_obs*2)=bundleData.obs(sort_idx*2);
        
    camMat=zeros(bundleData.n_cams, bundleData.n_cams);        
    for i=1:bundleData.n_pts
        n=size(pointTracks{i}, 1);
        for j=1:n
            for k=j+1:n
                row=pointTracks{i}(j);
                col=pointTracks{i}(k);
                if row>col
                    col=row;
                    row=pointTracks{i}(k);
                end
                camMat(row, col)=camMat(row, col)+1;
            end
        end
    end
    camMat=camMat+camMat';               
    fclose(fid);        
end