function normrs=read_normrs(filename, LM_iter)
    normrs=cell(LM_iter, 1);
    fid=fopen(filename, 'r');
    for i=1:LM_iter    
        numPCGIter=fscanf(fid, '%d', 1);
        normr=zeros(numPCGIter, 1);
        for j=1:numPCGIter
            input=fscanf(fid, '%f', 1);
            normr(j)=input;
        end
        normrs{i}=normr;
    end
    fclose(fid);
end