%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,obj_val_vec,iter_line_search,alpha_vec]=...
spectral_update_func(fun,x_k_cell,alpha_k,fun_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

    %%% Additional parameters used in globalization steps
    rho=0.8;
    M=10;gamma=10^(-4);

    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+alpha_k{1,i}.*fun_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
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

            LHS=obj_val_vec;%1*n_var
           
            if iter_line_search==1
                for i=1:n_var
                    d_k_cell{1,i}=x_k_plus_1_cell{1,i}-x_k_cell{1,i};
                end
            end

            for i=1:n_var 
                obj_val_PAST_MAX_i=max(obj_val_table(max(1,k+1-M+1):k+1,i));

                if spec.minimization_spec(1,i)==0         
                    eta_k_i=sqrt(obj_val_table(1,i))/((1+k)^2);
                    RHS_i=obj_val_PAST_MAX_i+eta_k_i-gamma*(max(alpha_k{i}(:)).^2)*(obj_val_table(k+1,i));%1*n_var

                else % minimization_spec==1
                    RHS_i=obj_val_PAST_MAX_i+gamma*(max(alpha_k{i}(:)))*sum(fun_k_cell{1,i}(:).*d_k_cell{1,i}(:));%1*n_var
                end %if spec.minimization_spec(1,i)==0
            end % for loop wrt i

            %%% If LHS<=RHS, exit the iteration
            continue_backtracing_dummy=(RHS_i-LHS<0);

            if sum(continue_backtracing_dummy(:))>0 & spec.positive_alpha_spec==1 
                for i=1:n_var
                    alpha_k{1,i}=alpha_k{1,i}.*(-rho);
                end

            elseif sum(continue_backtracing_dummy(:))>0 & (mod(iter_line_search,2)==0)  & spec.positive_alpha_spec==0
               for i=1:n_var
                  alpha_k{1,i}=alpha_k{1,i}.*(-rho);
              end

            elseif sum(continue_backtracing_dummy(:))==0 & mod(iter_line_search,2)==1 & spec.positive_alpha_spec==0
                for i=1:n_var
                   alpha_k{1,i}=alpha_k{1,i}*(-1);
                end
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

    for i=1:n_var
        alpha_vec(1,i)=max(alpha_k{i}(:)); 
    end % for loop wrt i

end
