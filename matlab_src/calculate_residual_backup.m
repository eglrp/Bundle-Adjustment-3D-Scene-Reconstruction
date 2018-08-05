% Outputs the function value with input vector state_x
function f=calculate_residual(state_x, bundleData) 
    f=zeros(bundleData.size_o*bundleData.n_obs, 1);
    for i=1:bundleData.n_obs
        camera=bundleData.cidx(i);
        point=bundleData.pidx(i);
        c=state_x(camera*bundleData.size_c+1:camera*bundleData.size_c+9);
        x=state_x(bundleData.n_cams*bundleData.size_c+point*bundleData.size_p+1:bundleData.n_cams*bundleData.size_c+point*bundleData.size_p+3);
        o=bundleData.obs((i-1)*bundleData.size_o+1:(i-1)*bundleData.size_o+2);
        fout=calculate_one_residual(c, x, o);                        
        f(i*2-1) = fout(1); 
        f(i*2) = fout(2);
    end
end
