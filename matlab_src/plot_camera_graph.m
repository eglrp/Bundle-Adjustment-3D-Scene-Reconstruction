function plot_camera_graph(A)
    radius=1.0;
    N=size(A, 1);
    xy=[radius*cos(0:2*pi/N:2*pi-2*pi/N)' radius*sin(0:2*pi/N:2*pi-2*pi/N)'];
    gplotd();
    
    function gplotd(varargin)
    %GPLOTD Plot a Directed Graph
    % GPLOTD(A,XY) Plots the Directed Graph represented by adjacency
    %   matrix A and points xy using the default style described below
    % GPLOTD(A,XY,PARAM1,VAL1,PARAM2,VAL2,...) Plots the Directed Graph
    %   using valid parameter name/value pairs
    %
    %   Inputs:
    %       A     - NxN adjacency matrix, where A(I,J) is nonzero
    %               if and only if there is an edge between points I and J
    %       xy    - Nx2 matrix of x/y coordinates
    %       ...   - Parameter name/value pairs that are consistent with
    %               valid PLOT parameters can also be specified
    %
    %   Default Plot Style Details:
    %       1. Undirected (2-way) edges are plotted in solid black lines
    %       2. Directed (1-way) edges are plotted in two styles
    %           a. If the edge connects a larger vertex ID with a smaller ID,
    %               the edge is plotted as a blue dashed line
    %           b. If the edge connects a smaller vertex ID with a larger ID,
    %               the edge is plotted as a red dotted line
    %       3. Any vertex that is connected to itself is plotted with a
    %           purple circle around it
    %
    %   Example:
    %       % plot a directed graph using default line styles
    %       n = 9; t = 2*pi/n*(0:n-1);
    %       A = round(rand(n));
    %       xy = 10*[cos(t); sin(t)]';
    %       gplotd(A,xy);
    %       for k = 1:n
    %           text(xy(k,1),xy(k,2),['  ' num2str(k)],...
    %               'Color','k','FontSize',12,'FontWeight','b')
    %       end
    % 
    %   Example:
    %       % plot a directed graph using plot name/value parameter pairs
    %       n = 9; t = 2*pi/n*(0:n-1);
    %       A = round(rand(n));
    %       xy = 10*[cos(t); sin(t)]';
    %       gplotd(A,xy,'LineWidth',2,'MarkerSize',8);
    %
    % See also: gplot, plot
    %
    % Author: Joseph Kirk
    % Email: jdkirk630@gmail.com
    % Release: 1.1
    % Release Date: 4/12/08

    % process inputs    
    [nr,nc] = size(A);
    [n,dim] = size(xy);
    if (~n) || (nr ~= n) || (nc ~= n) || (dim < 2)
        eval(['help ' mfilename]);
        error('Invalid input. See help notes above.');
    end
    params = struct();
    for var = 1:2:length(varargin)-1
        params.(varargin{var}) = varargin{var+1};
    end

    % parse the adjacency matrix
    A = double(logical(A));
    iA = diag(diag(A));         % self-connecting edges
    dA = A.*(1-A');             % directed edges (1-way)
    dAu = triu(dA,1);           % directed edges (1-way, sm2lg)
    dAl = tril(dA,-1);          % directed edges (1-way, lg2sm)
    uA = A-dA-iA;               % undirected edges (2-way)

    % make NaN-separated XY vectors
    [ix,iy] = makeXY(iA,xy);
    [ux,uy] = makeXY(tril(uA,0),xy);
    [dxu,dyu] = makeXY(dAu,xy);
    [dxl,dyl] = makeXY(dAl,xy);

    % plot the graph
    plot(ux,uy,'k-',params)
    hold on
    plot(dxu,dyu,'r:',params)
    plot(dxl,dyl,'b--',params)
    plot(ix,iy,'o','Color',[.6 0 1],params)
    plot(xy(:,1),xy(:,2),'k.')
    hold off

        function [x,y] = makeXY(A,xy)
            if any(A(:))
                [J,I] = find(A');
                m = length(I);
                xmat = [xy(I,1) xy(J,1) NaN(m,1)]';
                ymat = [xy(I,2) xy(J,2) NaN(m,1)]';
                x = xmat(:);
                y = ymat(:);
            else
                x = NaN;
                y = NaN;
            end
        end
    end
end