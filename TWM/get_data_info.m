function [data_max, data_min, exp_val, percentiles] = get_data_info(M, p);

% M = [5, 8, 1, 4; 2, 7, 3, 6];

if (isempty(p))
    p = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95];
end


num_elems = numel(M);
M_vect = reshape(M, [num_elems, 1]);

exp_val = mean(M_vect);
std_dev = std(M_vect);

z_vals = norminv(p);


data_max = max(M_vect);
data_min = min(M_vect);
percentiles = exp_val + z_vals*std_dev;



% perc = M_vect(p*num_elems)