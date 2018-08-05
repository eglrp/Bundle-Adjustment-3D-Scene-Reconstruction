% Point clouds center at (0, 0, 0), lies parallel to the z=0 plane, has a
% radius of 1. Number of layers indicates how many layers of such xy-plane
% rings are. Default is 1, in which case only one ring with z=0. Layer
% seperation indicates the z-distance between different layers. Layers stack upwards. Number of
% points indicate how many points are in each layer

% All cameras lie on the z=0 plane, sitting at radius r (r>1). Number of
% cameras indicate how many cameras are there sitting in the radius r ring.
% Assume all cameras have focal length f=100, and radial distortion
% parameters k1=k2=0.
function generate_synthetic_data(opts)
    num_cams=16;
    if isfield(opts, 'num_cams')
        num_cams=opts.num_cams;
    end
    cam_radius=10;
    if isfield(opts, 'cam_radius')
        cam_radius=opts.cam_radius;
    end
    num_pts_per_layer=48;
    if isfield(opts, 'num_pts_per_layer')
        num_pts_per_layer=opts.num_pts_per_layer;
    end
    num_layers=1;
    if isfield(opts, 'num_layers')
        num_layers=opts.num_layers;
    end
    layer_seperation=1;
    if isfield(opts, 'layer_seperation')
        layer_seperation=opts.layer_seperation;
    end
    num_pts_seen_per_cam_per_layer=6;
    if isfield(opts, 'num_pts_seen_per_cam_per_layer')
        num_pts_seen_per_cam_per_layer=opts.num_pts_seen_per_cam_per_layer;
    end
    output_filename='~/synthetic_data.txt';
    if isfield(opts, 'output_filename')
        output_filename=opts.output_filename;
    end
    copyCam0=false;
    if isfield(opts, 'copyCam0')
        copyCam0=true;
    end
    
    num_pts=num_pts_per_layer*num_layers;
    pt_positions=zeros(num_pts, 3);
    for i=1:num_layers
        pts_z=(i-1)*layer_seperation;
        offset=(i-1)*num_pts_per_layer;
        for j=1:num_pts_per_layer
            theta=(2*pi*j)/num_pts_per_layer;
            pt_positions(offset+j, 1)=cos(theta);
            pt_positions(offset+j, 2)=sin(theta);
            pt_positions(offset+j, 3)=pts_z+1;
        end
    end
        
    cam_positions=zeros(num_cams, 3);
    cam_orientations=zeros(num_cams, 3);
    cam_t=zeros(num_cams, 3);
    
    num_obs=num_cams*num_pts_seen_per_cam_per_layer*num_layers;
    obs_list=zeros(num_obs, 4);
    cp_cam_obs_list=[];
    pts_seen_by_cam1=[];
    obs_list_counter=1;
    for i=1:num_cams
        theta=(2*pi*i)/num_cams;
        cam_positions(i, 1)=cos(theta)*cam_radius;
        cam_positions(i, 2)=sin(theta)*cam_radius;        
        cam_positions(i, 3)=0.5;
   %     cam_positions(i, :)=cam_positions(i, :)+cam_positions(i, :)/10000.*[rand_func() rand_func() rand_func()];
        % v is the starting vector. Want to rotation v around k to [0 0
        % -1] for pi/2 radians
        v=[-cam_positions(i, 1) -cam_positions(i, 2) 0]/norm(cam_positions(i, :), 2);
        % k is axis of rotation
        k=cross(v, [0 0 -1]);
        cam_orientations(i, :)=k*pi/2;
        
        R=rodrigues(cam_orientations(i, :));
        t=-R*cam_positions(i, :)';
        cam_t(i, :)=t'; 
        
        dists=zeros(num_pts_per_layer, 1);
        for j=1:num_pts_per_layer
            dists(j)=norm(pt_positions(j, :)-cam_positions(i, :), 2);            
        end
        [~, idx]=sort(dists);
        for j=1:num_pts_seen_per_cam_per_layer
            for k=1:num_layers
                pt_idx=(k-1)*num_pts_per_layer+idx(j);
                obs_list(obs_list_counter, 1)=i;
                obs_list(obs_list_counter, 2)=pt_idx;
                p=calculate_projected_coordinates(pt_positions(pt_idx, :)', R, t, 100, 0, 0);
                obs_list(obs_list_counter, 3)=p(1);
                obs_list(obs_list_counter, 4)=p(2);
                obs_list_counter=obs_list_counter+1;
                
                if i==num_cams && sum(pts_seen_by_cam1==pt_idx)==0
                    R=rodrigues(cam_orientations(1, :));
                    t=cam_t(1, :)';
                    p=calculate_projected_coordinates(pt_positions(pt_idx, :)', R, t, 100, 0, 0);
                    cp_cam_obs_list=[cp_cam_obs_list; num_cams+1 pt_idx p(1) p(2)];                    
                elseif i==1
                    pts_seen_by_cam1=[pts_seen_by_cam1; pt_idx];
                end
            end
        end
    end
    
    fid=fopen(output_filename, 'w');
    if copyCam0
        fprintf(fid, '%d %d %d\n', num_cams+1, num_pts, num_obs+size(cp_cam_obs_list, 1));
    else
        fprintf(fid, '%d %d %d\n', num_cams, num_pts, num_obs);
    end
    for i=1:num_obs
        fprintf(fid, '%d %d %e %e\n', obs_list(i, 1)-1, obs_list(i, 2)-1, obs_list(i, 3), obs_list(i, 4));
    end
    if copyCam0
        for i=1:size(cp_cam_obs_list, 1)
            fprintf(fid, '%d %d %e %e\n', cp_cam_obs_list(i, 1)-1, cp_cam_obs_list(i, 2)-1, cp_cam_obs_list(i, 3)+rand_func()*obs_list(i, 3)/20, cp_cam_obs_list(i, 4)+rand_func()*obs_list(i, 4)/20);
        end
    end
    cp_cam_params=zeros(9, 1);
    noise_factor=999999999999999;
    for i=1:num_cams
        for j=1:3
            orient=cam_orientations(i, j)+rand_func()*cam_orientations(i, j)/noise_factor;
            fprintf(fid, '%10.13e\n', orient);
            if i==1
                cp_cam_params(j)=orient;
            end
        end
        for j=1:3
            t=cam_t(i, j)+rand_func()*cam_t(i, j)/noise_factor;
            if i==1
                cp_cam_params(j+3)=t;
            end
            fprintf(fid, '%10.13e\n', t);
        end
        focal=100;
        k1=0;
        k2=0;
        if i==1
            cp_cam_params(7)=focal;
            cp_cam_params(8)=k1;
            cp_cam_params(9)=k2;
        end
        fprintf(fid, '%10.13e\n', focal);
        fprintf(fid, '%10.13e\n', k1);
        fprintf(fid, '%10.13e\n', k2);        
    end
    if copyCam0
        for i=1:9
            fprintf(fid, '%10.13e\n', cp_cam_params(i));
        end
    end
    for i=1:num_pts
        for j=1:3
            fprintf(fid, '%10.13e\n', pt_positions(i, j)+rand_func()*pt_positions(i, j)/noise_factor);
        end
    end
    fclose(fid);
end

function out=rand_func()
    out=rand(1)-0.5;
end

function p2=calculate_projected_coordinates(X, R, t, f, k1, k2)
    P=R*X+t;
    p=-P/P(3);
    rp=1.0+k1*norm(p, 2)^2+k2*norm(p, 2)^4;
    p2=f*rp*p;
end