function mat = rmBDVal( mat )
%RMBDVAL Summary of this function goes here
%   Detailed explanation goes here
    s = size(mat);
    nrow = s(1) - 1; ncol = s(2) - 1;
    
    mat = mat(2:nrow, :);
    mat = mat(:, 2:ncol);
end

