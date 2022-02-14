function [H13]=significant_height(vectorh)
    n = length(vectorh);
    h = sort(vectorh,"ascend");
    H = 0;
    for i = round(2/3*n+1):n
        H = H + h(i);
                
    end
    
    H13 = 1/(1/3*n)*H;
    
    