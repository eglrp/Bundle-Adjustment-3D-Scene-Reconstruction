function [output, x]=read_bundle_output(filename)
    fid=fopen(filename, 'r');
    % Col 1: f value   Col 2: g value   Col 3: h   Col 4: rho    Col
    % 5: mu    Col 6: PCG flag   Col 7: PCG residual   Col 8: PCG
    % iteration
    output=[];
    line=fgets(fid);    %#ok<NASGU>
    while true
        line=fgets(fid);
        data=textscan(line, 'pcg flag: %d pcg residual: %f pcg iteration: %d');
        if size(data{1}, 1)==0
            data=textscan(line, '%d: f:%f g:%f h:%f rho:%f mu:%f', 1);
            if size(data{5}, 1)==0
                data=textscan(line, '%d: f:%f g:%f mu:%f', 1);
                if size(data{3}, 1)~=0
                    continue;
                else                    
                    break;
                end
            end
            output1=[data{2} data{3} data{4} data{5} data{6} 0 0 0];
            output=[output; output1];
        else
            line=fgets(fid);
            data=textscan(line, '%d: f:%f g:%f h:%f rho:%f mu:%f', 1);
            if size(data{6}, 1)==0
                disp('Corrupted data format when reading PCG\n');
                return;
            end
            output1=[data{2} data{3} data{4} data{5} data{6} 0 0 0]; %#ok<*AGROW>
            output=[output; output1];
        end        
    end
    %{
    for i=1:iteration   
        if isPCG
            line=fgets(fid);
            data=textscan(line, 'pcg flag: %d pcg residual: %f pcg iteration: %d');
            output(i, 6)=data{1};
            output(i, 7)=data{2};
            output(i, 8)=data{3};        
            line=fgets(fid);
            data=textscan(line, '%d: f:%f g:%f h:%f rho:%f mu:%f', 1);
            output(i, 1)=data{2};
            output(i, 2)=data{3};
            output(i, 3)=data{4};
            output(i, 4)=data{5};
            output(i, 5)=data{6};
        else
            line=fgets(fid);
            data=textscan(line, '%d: f:%f g:%f h:%f rho:%f mu:%f', 1);
            output(i, 1)=data{2};
            output(i, 2)=data{3};
            output(i, 3)=data{4};
            output(i, 4)=data{5};
            output(i, 5)=data{6};
        end
    end
    line=fgets(fid);
    %}
    x=fscanf(fid, '%f');
end