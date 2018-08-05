function plot_camera_evolution(cams, n_cams, iter_start, iter_end, varargin)
    cams=cams(n_cams*9*(iter_start-1)+1:n_cams*9*iter_end);
    n=iter_end-iter_start+1;
    
    cam_pos=zeros(n*n_cams, 3);
    cam_dir=zeros(n*n_cams, 3);    
    for i=1:n
        x=cams((i-1)*n_cams*9+1:i*n_cams*9);
        for j=1:n_cams
            rotations=x(j*9-8:j*9-6);
            R=rodrigues(rotations);
            translations=x(j*9-5:j*9-3);
            start_pos=-R'*translations;
            end_pos=R'*([0;0;1]-translations);
            cam_pos((i-1)*n_cams+j, :)=start_pos';
            cam_dir((i-1)*n_cams+j, :)=(end_pos'-start_pos')/norm(end_pos'-start_pos');
        end
    end
    
    cam_pos_sum=zeros(n_cams, 3);
    for i=1:n_cams
        cam_pos_sum(i, :)=sum(cam_pos(i:n_cams:size(cam_pos, 1), :));
    end
    cam_pos_sum_ave=sum(cam_pos_sum)/size(cam_pos_sum, 1);
    
    drop_cams=[];
    % Drop some extreme cameras
    for i=1:5
        [~, max_idx]=max(cam_pos_sum);
        [~, min_idx]=min(cam_pos_sum);
        next_idx=[max_idx';min_idx'];
        drop_cams=[drop_cams; next_idx];
        cam_pos_sum([max_idx'; min_idx'], :)=ones(size(next_idx, 1), 1)*cam_pos_sum_ave;
    end
    
    drop_indices=[];
    for i=1:0
        offset=n_cams*(i-1);
        for j=1:size(drop_cams, 1)
            drop_indices=[drop_indices; offset+drop_cams(j)];
        end
    end
    cam_pos(drop_indices, :)=[];
    cam_dir(drop_indices, :)=[];
    n_cams=size(cam_pos, 1)/n;
    
    maxscale=max(cam_pos);
    minscale=min(cam_pos); 
    
    %{
    scale=(abs(maxscale)+abs(minscale))/2;
    cam_pos=((cam_pos-ones(size(cam_pos, 1), 1)*(maxscale+minscale)/2));
    cam_pos(:, 1)=cam_pos(:, 1)/scale(1);
    cam_pos(:, 2)=cam_pos(:, 2)/scale(2);
    cam_pos(:, 3)=cam_pos(:, 3)/scale(3);
        axis([-1.5 1.5 -1.5 1.5]);
    %}
    
    
    camera_pos=(maxscale+minscale)/2;
    camera_pos(3)=maxscale(3)*1.1;    
    camera_target=(maxscale+minscale)/2;
    
    colors=[abs(sin(linspace(0, 1, n_cams)*pi*30))' linspace(1, 0.5, n_cams)' abs(cos(linspace(0, 1, n_cams)*pi*10))'];    
    for i=1:n
        cam_pos2=cam_pos((i-1)*n_cams+1:i*n_cams, :);
        cam_dir2=cam_dir((i-1)*n_cams+1:i*n_cams, :);
        clf
        hold on;        
        campos('manual');
        campos(camera_pos);
        camtarget('manual');
        axis([minscale(1) maxscale(1) minscale(2) maxscale(2)]);
        length=sqrt((maxscale(1)-minscale(1))^2+(maxscale(2)-minscale(2))^2);
        alpha=-length/33;
        camtarget(camera_target);
   %     arrow('Start', cam_pos2, 'Stop', cam_dir2*0.07+cam_pos2, 'TipAngle', 3);  
        for j=1:size(cam_pos2)
            plot3(cam_pos2(j, 1), cam_pos2(j, 2), cam_pos2(j, 3), '.', 'Color', colors(j, :));  
            plot3([cam_pos2(j, 1); cam_pos2(j, 1)+alpha*cam_dir(j, 1)], [cam_pos2(j, 2); cam_pos2(j, 2)+alpha*cam_dir(j, 2)], [cam_pos2(j, 3); cam_pos2(j, 3)+alpha*cam_dir(j, 3)], '-k');            
        end        
        hold off;         
        M(i)=getframe(gcf);   
    end
    
    if nargin>=5
        movie2avi(M, varargin{1}, 'fps', 1);
    end
end