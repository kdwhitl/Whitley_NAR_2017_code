function mat_out = cutmat(mat,column,range)
%cut a matrix by the values in one column specified by range ([min max])
%100426, mjc

[~, c] = size(mat);

mskcol = mat(:,column) >= range(1) & mat(:,column) <= range(2);

mskmat = repmat(mskcol,1,c);

mat_out = reshape(mat(mskmat),[],c);

end

