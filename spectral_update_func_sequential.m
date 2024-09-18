%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
obj_val_vec,step_size]=...
spectral_update_func_sequential(fun,x_k_cell,alpha_k,d_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

global alpha_k_original

    step_size=ones(1,n_var);
    alpha_k_original=alpha_k;
    obj_val_vec=[];
    merit_func=spec.merit_func;

    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
        %%%%%%%%%%%%%%%%%%
        if spec.merit_func_spec==1 & iter_line_search==2%%%%%%
            d_k_cell{1,1}=d_k_cell{1,1}./alpha_k{1,1};
            alpha_k{1,1}=1;%%%%%%
        end%%%%%
        %%%%%%%%%%%%%%%%%%

       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+step_size(i)*d_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
       end

        if spec.line_search_spec==1

                merit_obj_k=obj_val_table(k+1,1);


                x_k_plus_1_cell_temp=x_k_plus_1_cell;
                x_k_plus_1_cell_temp{1,2}=x_k_cell{1,2};
                x_k_plus_1_cell_temp{1,3}=x_k_cell{1,3};
                
                if isempty(spec.other_input_merit_func)==1
                    merit_obj_k_plus_1=merit_func(x_k_plus_1_cell_temp);
                else
                    merit_obj_k_plus_1=merit_func(x_k_plus_1_cell_temp,spec.other_input_merit_func);
                end

                obj_val_vec=merit_obj_k_plus_1;

                merit_obj_k_plus_1=merit_obj_k_plus_1;

                %%% Simple
                %eta_k=0.0001/((1+k)^2);
                eta_k=0.1/((1+k)^2);
                %eta_k=10/((1+k)^2);%%%%
                %eta_k=1000/((1+k)^2);%%%%
                
                LHS=merit_obj_k_plus_1;
                RHS=merit_obj_k+eta_k;
                %[LHS,RHS]
                continue_backtracking_dummy=(LHS>=RHS);

            rho=spec.rho;


            if continue_backtracking_dummy==1 & spec.positive_alpha_spec==1

                %%%%%%%%
                if spec.merit_func_spec==1 & iter_line_search==2
                    step_size(1)=0;%%%%%%
                    obj_val_vec=[];
                end
                %%%%%%%%%%%
                
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
        %alpha_k{1,1}
        %alpha_vec(1,i)=max(alpha_k_original{i}(:)); 
        alpha_vec(1,i)=max(alpha_k{1,i}(:)); 
    end


    if step_size(1)==0% If no update of parameters, update x and dx_dparam using the parameters
        x_k_plus_1_cell{1,1}=x_k_cell{1,1};
        x_k_plus_1_cell{1,2}=x_k_cell{1,2};
        x_k_plus_1_cell{1,3}=x_k_cell{1,3};
        other_input_cell{end-1}=2;%%% Compute x_updated, dx_dparam_updated, using param_updated and x_initial,dx_dparam_initial
        [x_k_plus_1_cell,~]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});
    end
  

        
        if spec.bound_spec==1
            x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
        end
        
        other_input_cell{end-1}=0;% Compute param_updated,x_updated,dx_dparam_updated
        [x_k_plus_2_cell,other_output_k_plus_1]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});
       


             for i=1:n_var
                 fun_k_plus_1_cell{1,i}=x_k_plus_2_cell{1,i}-x_k_plus_1_cell{1,i};
             end


        %%% DIST: sup norm of F(x)=x-Phi(x). 
        DIST_vec=ones(1,n_var);
        for i=1:n_var
            DIST_vec(1,i)=norm_func(fun_k_plus_1_cell{1,i}(:),x_k_plus_1_cell{1,i}(:),spec.norm_spec(i));
        end

        if isempty(obj_val_vec)==1
            
            if isempty(spec.other_input_merit_func)==1
                obj_val_vec=merit_func(x_k_plus_1_cell);
            else
                obj_val_vec=merit_func(x_k_plus_1_cell,spec.other_input_merit_func);
            end
        end

%%%%%%%%%%%%%%%%%

end
