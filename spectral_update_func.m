%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
obj_val_vec,step_size]=...
spectral_update_func(fun,x_k_cell,alpha_k,fun_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

    global step_size

    
    for i=1:n_var      
        d_k_cell{1,i}=alpha_k{1,i}.*fun_k_cell{1,i};
    end

    step_size=ones(1,n_var);


    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+step_size(i)*d_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
       end

     if spec.merit_func_spec==0
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

    end % merit_func_spec==0

        if spec.line_search_spec==1

           if spec.merit_func_spec==0

              %%% Implicitly assume common step size for all the variables %%%%
              continue_backtracking_dummy=line_search_terminate_func(...
                       obj_val_vec,obj_val_table,n_var,k,step_size(1),spec);

            else% merit_func_spec==1
                merit_func=spec.merit_func;
                merit_obj_k=obj_val_table(k+1,1);
                if isempty(spec.other_input_merit_func)==1
                    merit_obj_k_plus_1=merit_func(x_k_plus_1_cell);
                else
                    merit_obj_k_plus_1=merit_func(x_k_plus_1_cell,spec.other_input_merit_func);
                end

                obj_val_vec=merit_obj_k_plus_1;

                merit_obj_k_plus_1=merit_obj_k_plus_1;
                %%% Simple
                eta_k=0/((1+k)^2);
                LHS=merit_obj_k_plus_1+eta_k;
                RHS=merit_obj_k;
                %%[LHS,RHS]
                continue_backtracking_dummy=(LHS>=RHS);
                
                if 1==0
                %%% Similar to Cruz et al. 2006
                

                M=spec.M;
                gamma=spec.gamma;
             
                obj_val_table_sum=sum(obj_val_table,2);
                LHS=sum(obj_val_vec);%1*n_var
               
                obj_val_PAST_MAX=max(obj_val_table_sum(max(1,k+1-M+1):k+1,:));
                eta_k=sqrt(obj_val_table_sum(1,:))/((1+k)^2);
                
                eta_k=1*(obj_val_table_sum(1,:))/((1+k)^2);%%%%
            
                d_norm_squared=0;
                for i=1:n_var
                    d_norm_squared=d_norm_squared+sum(d_k_cell{1}(:).^2);
                end

                RHS=obj_val_PAST_MAX+eta_k-(gamma*step_size(1)^2)*(d_norm_squared);


                continue_backtracking_dummy=(RHS-LHS<0);
                end

           end

            rho=spec.rho;



            if continue_backtracking_dummy==1 & spec.positive_alpha_spec==1
                %step_size=step_size.*rho; 
                %step_size=step_size.*(-rho); 
                
                %%%%%
                step_size(1)=0;%%%%%%
                
                
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
        alpha_vec(1,i)=max(alpha_k{i}(:)); 
    end

%%%%%%%%%%%%%
     if spec.merit_func_spec==1
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

    end
%%%%%%%%%%%%%%%%%

end
