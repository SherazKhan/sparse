function log_det = find_log_det(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Camilo Lamus (lamus@mit.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_det = log(A(1,1));
n = size(A);
n = n(1,1);
for i = 1:(n-1)
    A = sweep(A);
    log_det = log_det + log(A(1,1));
    log(A(1,1)); 
end

end