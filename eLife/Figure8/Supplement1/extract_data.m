%% Extract data from files
%
% base_dir - base data direcotry (e.g. Domain5940_142). Ot assumest that
% under this directory, there a directory called 'tensor_store' with all .m
% files.
%

function [M0, M, V, U] = extract_data(base_dir)

    M0 = [];
    M = [];
    V = [];
    U = [];
    files = dir(base_dir + filesep + "tensor_store" + filesep + "*.mat");
    for i = 1:size(files,1)
        tensor_store = load(base_dir + filesep + "tensor_store" + filesep + files(i).name, "-mat", "tensor_store");
        ts = struct2table(tensor_store.tensor_store);
        M0 = [M0; ts.M0(:,1:4)];
        M = [M; ts.M(:,1:4)./ts.M(:,5)];
        V = [V; ts.V(:,1:4)];
        U = [U; ts.dUdt(:,1:4)];
    end

end