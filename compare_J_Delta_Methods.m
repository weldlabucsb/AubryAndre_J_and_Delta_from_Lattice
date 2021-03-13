function outTable = compare_J_Delta_Methods(s1,s2,lambda1,lambda2)
tic
    if size(s1,2)==1
        s1 = transpose(s1);
    end
    if size(s2,2)==1
        s2 = transpose(s2);
    end
    if size(lambda1,2)==1
        lambda1 = transpose(lambda1);
    end
    if size(lambda2,2)==1
        lambda2 = transpose(lambda2);
    end
    
    s1s = [];
    s2s = [];
    la1s = [];
    la2s = [];
    
    J_gauss = [];
    J_num_approx = [];
    J_num = [];
    J_piecewise = [];
    
    De_gauss = [];
    De_num = [];
    De_piecewise = [];
    
    for s1_i = s1
        for s2_i = s2
            for la1_i = lambda1
                for la2_i = lambda2
                    s1s = [s1s;s1_i];
                    s2s = [s2s;s2_i];
                    la1s = [la1s;la1_i];
                    la2s = [la2s;la2_i];
                    
                    [J_gauss_i,De_gauss_i] = J_Delta_Gaussian(s1_i,s2_i,la1_i,la2_i);
                    J_gauss = [J_gauss; J_gauss_i];
                    De_gauss = [De_gauss; De_gauss_i];
                    
                    J_num_approx_i = J_Numeric_Approx(s1_i);
                    J_num_approx = [J_num_approx;J_num_approx_i];
                    
                    [J_num_i,De_num_i] = J_Delta_Numeric(s1_i,s2_i,la1_i,la2_i,0);
                    J_num = [J_num;J_num_i];
                    De_num = [De_num;De_num_i];
                    
                    [J_piecewise_i, De_piecewise_i] = J_Delta_PiecewiseFit(s1_i,s2_i);
                    J_piecewise = [J_piecewise;J_piecewise_i];
                    De_piecewise = [De_piecewise;De_piecewise_i];

                end
            end
        end
    end
    
    tabEntries = {'s1','s2','lambda1','lambda2','J Gauss W Approx', 'J Approx Formula','J Numeric', 'J Piecewise Fit', 'Delta Gauss W Approx','Delta Numeric','Delta Piecewise Fit'};
    outTable = table(s1s,s2s,la1s,la2s,J_gauss,J_num_approx,J_num,J_piecewise,De_gauss,De_num,De_piecewise,'VariableNames',tabEntries);
toc
end

