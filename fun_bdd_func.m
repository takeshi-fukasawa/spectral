function [x_cell,fun_cell,other_output]=...
    fun_bdd_func(fun,x_cell,other_input_cell,x_max_cell,x_min_cell,n_var)
 %% Restrict the range of x
 
   for i=1:n_var
        if isempty(x_max_cell{1,i})==0
             x_cell{1,i}=min(x_cell{1,i},x_max_cell{1,i});
            ub_dummy_cell{1,i}=(x_cell{1,i}==x_max_cell{1,i});% upper bound dummy
        end
        if isempty(x_min_cell{1,i})==0
             x_cell{1,i}=max(x_cell{1,i},x_min_cell{1,i});
             lb_dummy_cell{1,i}=(x_cell{1,i}==x_min_cell{1,i});% lower bound dummy
       end
   end% for loop

 %% At boundary => Set fun val to 0
   [fun_cell,other_output]=...
           fun(x_cell{:},other_input_cell{:});

   for i=1:n_var
        if isempty(x_max_cell{1,i})==0
            %%% Upper bound => Set the func value to zero
             fun_cell{1,i}=fun_cell{1,i}.*(1-ub_dummy_cell{1,i});
        end
        if isempty(x_min_cell{1,i})==0
            %%% Lower bound => Set the func value to zero
             fun_cell{1,i}=fun_cell{1,i}.*(1-lb_dummy_cell{1,i});
        end
   end% for loop

end
