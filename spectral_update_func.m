%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
obj_val_vec,step_size]=...
spectral_update_func(fun,x_k_cell,fun_k_cell,alpha_k,d_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

    step_size=ones(1,n_var);
    
    %%% Update variables in the spectral algorithm %%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+step_size(i)*d_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
       end

        [fun_k_plus_1_cell,other_output_k_plus_1]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});
       

        if spec.fixed_point_iter_spec==1
             for i=1:n_var
                 fun_k_plus_1_cell{1,i}=fun_k_plus_1_cell{1,i}-x_k_plus_1_cell{1,i};
             end
        end
 
        %%% DIST: sup norm of F(x)=x-Phi(x). 
        DIST_vec=ones(1,n_var);
        for i=1:n_var
            DIST_vec(1,i)=norm_func(fun_k_plus_1_cell{1,i}(:),x_k_plus_1_cell{1,i}(:),spec.norm_spec(i));
        end
        obj_val_vec=DIST_vec.^2;


        if spec.line_search_spec==1

            %%% Implicitly assume common step size for all the variables %%%%
            continue_backtracking_dummy=line_search_terminate_func(...
                     obj_val_vec,obj_val_table,n_var,k,step_size(1),spec);


            rho=spec.rho;


            if continue_backtracking_dummy==1 & spec.positive_alpha_spec==1
                step_size=step_size.*rho; 
                
            elseif continue_backtracking_dummy==1 & (mod(iter_line_search,2)==0)  & spec.positive_alpha_spec==0
                step_size=step_size.*(-rho);

            elseif continue_backtracking_dummy==1 & mod(iter_line_search,2)==1 & spec.positive_alpha_spec==0
                step_size=step_size.*(-1);
                 
            else
                break;
            end

        end% line search spec==1
    end % end iter_line_search=1:ITER_MAX_LINE_SEARCH loop

    for i=1:n_var
        alpha_vec(1,i)=max(alpha_k{1,i}(:));     
    end


end
