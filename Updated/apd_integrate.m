function out = apd_integrate(mat,factor)
%take a vector and integrate it (like increasing APD integration time)
%090113 mjc

N = length(mat);
Ndesired = floor(N/factor)*factor;

if Ndesired < N
    mat = mat(1:Ndesired); 
end

temp = reshape(mat,factor,floor(N/factor));

if factor > 1
    out = sum(temp);
else
    out = temp;

end