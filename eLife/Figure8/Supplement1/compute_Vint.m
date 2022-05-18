%% Compute integrated V
%
% M0 - vector with componets of the M0 tensor. Rows are different time points,
% colums are M0_xx, M0_xy, M0_yx, M0_yy (sometimes ignored)
% M - vector with componets of the M tensor. Rows are different time points,
% colums are M_xx, M_xy, M_yx, M_yy 
% V - vector with componets of the V tensor. Rows are different time points,
% colums are V_xx, V_xy, V_yx, V_yy 
% M0_type - string
%    "experiment" - use M0 read from experimental data
%    "Msym" -  use M at the first time point; compute its eigenvalues (ev1
%              and ev2) and define M0 = diag([0.5*(ev1 + ev2), 0.5*(ev1 +
%              ev2)]); M0 passed to the function is ignored
%    "unit" -  set M0 = daig([1 1]); M0 passed to the function is ignored
%    "M0sym" - same as "Msym" but use M0 at the first time point;
%              M0(2:end,:) is ignored.
%    otherwise - M0 is set to M0 at the first time step, i.e. M0(1,:); 
%                M0(2:end,:) is ignored. 
%
function Vint = compute_Vint(M0, M, V, MO_type)

    use_exp_M0 = false;
    switch MO_type
        case "experiment"
            use_exp_M0 = true;
        case "Msym"
            ev = eig(reshape(M(1,:), [2,2]));
            MM0 = diag([0.5*mean(ev), 0.5*mean(ev)]);
        case "unit"
            MM0 = diag([1, 1]);
        case "M0sym"
            ev = eig(reshape(M0(1,:), [2,2]));
            MM0 = diag([0.5*mean(ev), 0.5*mean(ev)]);
        otherwise
            MM0 = reshape(M0(1,:), [2,2]);
    end
    

    U = zeros(size(V));
    for i = 1:size(V,1)
        if use_exp_M0
            MM0 = reshape(M0(i,:), [2,2]);
        end
        MM  = reshape(M(i,:), [2,2]);
        %U(i,:) = reshape(compute_U(MM0, MM), 1, []);
    end

    VV = zeros(size(V));
    %meanU = mean(U);
    %meanU = reshape(meanU, [2,2]);
    %[UU, D] = eig(meanU, 'vector');
    %[D, ind] = sort(D, 'descend');
    %UU = UU(:, ind);
    UU = [0 1; -1 0];
    for i = 1:size(V,1)
        VVtemp = reshape(V(i,:), [2,2]);
        VVtemp = UU' * VVtemp * UU;
        VV(i,:) = reshape(VVtemp, 1, []);
    end
    Vint = cumsum(VV);
end