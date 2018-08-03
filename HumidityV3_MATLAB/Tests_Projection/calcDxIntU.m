function res_vec = calcDxIntU(intU_vec, Dx)
%CALCDXINTU Calculate the x-derivative of the int_pA^pB u dp
%   INPUTS:: intU_vec: A vector of size Nx. Its ith conponent is the value int_pA^pB u(x_i,p)dp
%   OUTPUT:: res_vec: A vector of size Nx. Its ith conponent is the value (d/dx)int_pA^pB u(x_i,p)dp

resLen = length(intU_vec) - 1;  % Length of the result vector
res_vec = zeros(resLen, 1);
for i = 1:resLen
    res_vec(i) = (intU_vec(i + 1) - intU_vec(i)) / Dx;
end

end