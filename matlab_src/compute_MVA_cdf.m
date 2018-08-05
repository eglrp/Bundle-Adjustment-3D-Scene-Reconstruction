function angles_cum=compute_MVA_cdf(angles, bins)
    angles_cum=zeros(size(bins));
    for i=1:size(angles, 1)
        if angles(i)>bins(50)
            continue;
        end
        indices=find(bins>angles(i));
        angles_cum(indices)=angles_cum(indices)+1;
    end
end

function angles=find_max_angles(G, state_x)    
    angles=zeros(size(G, 2), 1);
    for i=1:size(G, 2)
        angle=find_max_angle(i, G, state_x, 450);
        angles(i)=angle;
    end
end

function maxangle=find_max_angle(pt_idx, G, state_x, n_cams)
    out=[];
    cams=find(G(:, pt_idx)>0);
    for i=1:size(cams, 1)
        cam=cams(i);
        cam_params=state_x(cam*9-8:cam*9);    
        rotations=cam_params(1:3);
        R=rodrigues(rotations);
        translations=state_x(4:6);
        cam_pos=-R'*translations;
        out=[out; cam_pos'];
    end
    point=state_x(n_cams*9+pt_idx*3-2:n_cams*9+pt_idx*3)';
    maxangle=0;
    for i=1:size(out, 1)
        for j=1:i-1
            v1=out(i, :)-point;
            v2=out(j, :)-point;
            angle=acos((v1*v2')/(norm(v1, 2)*norm(v2, 2)));
            if angle>maxangle
                maxangle=angle;
            end
        end
    end
end