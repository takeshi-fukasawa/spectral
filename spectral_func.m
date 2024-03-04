function [x_sol_cell,other_output_k,DIST_table,iter_info,fun_k_cell]=...
    spectral_func(fun,spec,...
    x_0_cell,varargin)

    n_var=size(x_0_cell,2);

    x_max_cell=cell(1,n_var);
    x_min_cell=cell(1,n_var);

   
    [x_sol_cell,other_output_k,DIST_table,iter_info,fun_k_cell]=...
        spectral_bound_func(fun,spec,...
        x_0_cell,x_max_cell,x_min_cell,varargin{:});

end
