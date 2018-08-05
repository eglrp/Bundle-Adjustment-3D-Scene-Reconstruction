function plot_reconstruction(x, n_cams, varargin)
    n_pts=(size(x, 1)-n_cams*9)/3;
    
    cam_pos=zeros(n_cams, 3);
    cam_dir=zeros(n_cams, 3);
    colors=[zeros(n_cams, 1) linspace(1, 0, n_cams)' linspace(0, 1, n_cams)'];
    for i=1:n_cams
        rotations=x(i*9-8:i*9-6);
        R=rodrigues(rotations);
        translations=x(i*9-5:i*9-3);
        start_pos=-R'*translations;
        end_pos=R'*([0;0;1]-translations);
        cam_pos(i, :)=start_pos';
        cam_dir(i, :)=(end_pos'-start_pos')/norm(end_pos'-start_pos');        
    end  
    
    arrow('Start', cam_pos, 'Stop', -cam_dir*1.2+cam_pos, 'TipAngle', 3);  
    for i=1:n_cams
        plot3(cam_pos(i, 1), cam_pos(i, 2), cam_pos(i, 3), '.', 'Color', colors(i, :));  
    end
    
    if n_pts>0
        point_pos=zeros(1, 3);
        offset=n_cams*9;
        for i=1:n_pts
            point_pos(i, :)=x(offset+i*3-2:offset+i*3)';
        end
        plot3(point_pos(:, 1), point_pos(:, 2), point_pos(:, 3), '.y');
    end
   
    % produce an off file to read into meshlab
    if nargin==3
        fid=fopen(varargin{1}, 'w');        
        
        fprintf(fid, 'COFF\n');        
        fprintf(fid, '%d %d %d\n', n_cams+n_pts, size(cam_pos, 1), 0);
        for i=1:n_cams
            fprintf(fid, '%d %d %d %d %d %d %d\n', cam_pos(i, 1), cam_pos(i, 2), cam_pos(i, 3), 0, 1, 0, 1);  
            fprintf(fid, '%d %d %d %d %d %d %d\n', cam_pos(i, 1)+cam_dir(i, 1)*0.07+0.01, cam_pos(i, 2)+cam_dir(i, 2)*0.07+0.01, cam_pos(i, 3)+cam_dir(i, 3)*0.07+0.01, 0, 0.3, 0, 1); 
            fprintf(fid, '%d %d %d %d %d %d %d\n', cam_pos(i, 1)+cam_dir(i, 1)*0.07-0.01, cam_pos(i, 2)+cam_dir(i, 2)*0.07-0.01, cam_pos(i, 3)+cam_dir(i, 3)*0.07-0.01, 0, 0.3, 0, 1); 
        end
        offset=n_cams*9;
        for i=1:n_pts
            fprintf(fid, '%d %d %d %d %d %d %d\n', x(offset+i*3-2), x(offset+i*3-1), x(offset+i*3), 1, 1, 1, 1);            
        end   
        for i=1:n_cams
            fprintf(fid, '%d %d %d %d %d %d %d %d\n', 3, i*3-3, i*3-2, i*3-1, 0, 0.5, 0, 0.7);
        end
        
        %{
        fprintf(fid, 'OFF\n');
        fprintf(fid, '%d %d %d\n', n_pts, 0, 0);
        offset=n_cams*9;
        for i=1:n_pts
            fprintf(fid, '%d %d %d\n', x(offset+i*3-2), x(offset+i*3-1), x(offset+i*3));            
        end
        %}
        fclose(fid);
    end
end    