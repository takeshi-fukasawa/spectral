%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,obj_val_vec,iter_line_search,step_size]=...
spectral_update_func(fun,x_k_cell,fun_k_cell,d_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

    %%% Additional parameters used in globalization steps
    rho=0.1;
    M=5;gamma=10^(-4);
    step_size=1;

    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+step_size*d_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1 & iter_line_search==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
         if spec.line_search_spec==1
        for i=1:n_var
           d_k_cell{1,i}=x_k_plus_1_cell{1,i}-x_k_cell{1,i};
       end
        end

       end

        if sum(spec.minimization_spec(:))==0 % Solve nonlinear eq
           [fun_k_plus_1_cell,other_output_k_plus_1]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});

            if spec.fixed_point_iter_spec==1
                for i=1:n_var
                    if spec.minimization_spec(1,i)==0
                        fun_k_plus_1_cell{1,i}=fun_k_plus_1_cell{1,i}-x_k_plus_1_cell{1,i};
                    end
                    
                end
            end

            %%% DIST: sup norm of F(x)=x-Phi(x). 
            DIST_vec=NaN(1,n_var);
            for i=1:n_var
                DIST_vec(1,i)=norm_func(fun_k_plus_1_cell{1,i}(:),x_k_plus_1_cell{1,i}(:),spec.norm_spec(i));
            end
            obj_val_vec=DIST_vec.^2;

        else % Minimization/ nonlinear eq spec

            [obj_fun_k_plus_1_cell]=...
              fun(x_k_plus_1_cell{:},other_input_cell{:});

            for i=1:n_var
                obj_val_vec(1,i)=obj_fun_k_plus_1_cell{1,i};
            end
        end%if else
       
   
        if spec.line_search_spec==1


continue_backtracking_dummy=line_search_terminate_func(obj_val_vec,obj_val_table,n_var,k,M,gamma,step_size,fun_k_cell,d_k_cell,spec);

            if sum(continue_backtracking_dummy(:))>0 & spec.positive_alpha_spec==1
                    step_size=step_size.*rho; 

            elseif sum(continue_backtracking_dummy(:))>0 & (mod(iter_line_search,2)==0)  & spec.positive_alpha_spec==0
                    step_size=step_size.*(-rho); 

            elseif sum(continue_backtracking_dummy(:))>0 & mod(iter_line_search,2)==1 & spec.positive_alpha_spec==0
                    step_size=step_size.*(-1); 
            else
                break;
            end

        end% line search spec==1
    end % end iter_line_search=1:ITER_MAX_LINE_SEARCH loop


    if sum(spec.minimization_spec(:))>0
        [obj_fun_k_plus_1_cell,fun_k_plus_1_cell,other_output_k_plus_1]=...
            fun(x_k_plus_1_cell{:},other_input_cell{:});

        for i=1:n_var
            obj_val_vec(1,i)=obj_fun_k_plus_1_cell{1,i};
            DIST_vec(1,i)=norm_func(fun_k_plus_1_cell{1,i}(:),x_k_plus_1_cell{1,i}(:),spec.norm_spec(i));
        end

    end

end
