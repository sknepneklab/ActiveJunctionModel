%% compute U tensor
%
%
%  source: eq (14) in F. Franer, et al., Eur. Phys. J. E 25, 349 (2008)
%  M  - a 2x2 matrix
function U = compute_U(M)

    meanM = reshape(mean(M), [2,2]);
    [UU, D] = eig(meanM, 'vector');
    [D, ind] = sort(D, 'descend');
    meanD = mean(D);
    M0 = diag([meanD, meanD]);

    U = zeros(size(M));
    for i = 1:size(M,1)
        MM = reshape(M(i,:), [2,2]);
        [UU, D] = eig(MM, 'vector');
        [D, ind] = sort(D, 'descend');
        MM = UU' * diag(D) * UU;
        U(i,:) = reshape(0.5*(logm(MM) - logm(M0)), 1, []);
    end
end