function [Hrms]=rms_height(vectorh)
    n = length(vectorh);
    h = sort(vectorh,"descend");
    
    Hrms = 0;
    for i = 1:n
        Hrms = Hrms + h(i).^2; 
    end

    Hrms = sqrt(1/n*Hrms);
