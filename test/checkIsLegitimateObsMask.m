function good=checkIsLegitimateObsMask(bundleData, camMat, obs_mask)
    good=true;
    pointTracks2=cell(bundleData.n_pts, 1);
    camCount=zeros(bundleData.n_cams, 1);
    for i=1:size(obs_mask, 1)
        obs_idx=obs_mask(i);
        cid=bundleData.cidx(obs_idx);
        pid=bundleData.pidx(obs_idx);
        pointTracks2{pid+1}=[pointTracks2{pid+1}; cid+1];
        camCount(cid+1)=camCount(cid+1)+1;
    end
        
    for i=1:bundleData.n_cams
        if camCount(i)==0
            disp('Some camera has been dropped out completely');
            good=false;
            return;
        end
    end
    for i=1:size(pointTracks2, 1)
        if size(pointTracks2{i}, 1)==0
            good=false;
            disp('Some point has no camera seeing it');
            return;
        end
    end
    
    camMat2=zeros(bundleData.n_cams, bundleData.n_cams);
    for i=1:bundleData.n_pts
        n=size(pointTracks2{i}, 1);
        for j=1:n
            for k=j+1:n
                row=pointTracks2{i}(j);
                col=pointTracks2{i}(k);
                if row>col
                    col=row;
                    row=pointTracks2{i}(k);
                end
                camMat2(row, col)=camMat2(row, col)+1;
            end
        end
    end
    camMat2=camMat2+camMat2';
    camMat=camMat+camMat';
    camMat=camMat-diag(diag(camMat));
    
    [row, col]=find(camMat~=0);
    [row2, col2]=find(camMat2~=0);
    M1=sparse(row, col, ones(size(row, 1), 1));
    M2=sparse(row2, col2, ones(size(row2, 1), 1));
    out=(find(M1-M2)<0);
    if sum(out)>0
        disp('Resulting graph is not a subset of input skeletal graph');
        good=false;
    end
end