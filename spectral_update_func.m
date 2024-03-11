%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec]=...
spectral_update_func(fun,x_k_cell,alpha_k,fun_k_cell,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,DIST_table)

    rho=0.8;
    M=10;gamma=10^(-4);

    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}-alpha_k{1,i}.*fun_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);
     end

        [fun_k_plus_1_cell,other_output_k_plus_1]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});
       

        if spec.fixed_point_iter_spec==1
             for i=1:n_var
                 fun_k_plus_1_cell{1,i}=x_k_plus_1_cell{1,i}-fun_k_plus_1_cell{1,i};
             end
        end
 
        %%% DIST: sup norm of F(x)=x-Phi(x). 
        DIST_vec=ones(1,n_var);
        for i=1:n_var
            DIST_vec(1,i)=max(abs(fun_k_plus_1_cell{1,i}),[],'all','omitnan');
        end
   
        if spec.line_search_spec==1

            id=1;
            if spec.iter_line_search==3
                alpha_k{1,id}=alpha_k{1,id}*1;
            end

            DIST_PAST_MAX_vec=max(DIST_table(max(1,k+1-M+1):k+1,:),[],1);

            eta_k=DIST_table(1,id)/((1+k)^2);
            LHS=DIST_vec(id);

            %%%%% alpha_x^2+alpha_y^2???
            RHS=DIST_PAST_MAX_vec(id)+eta_k-gamma*(sum(alpha_k{1,id}(:).^2))*DIST_table(k+1,id);

            if LHS>RHS & mod(iter_line_search,2)==0
                alpha_k{1,id}=alpha_k{1,id}.*(-rho);
            elseif LHS>RHS & mod(iter_line_search,2)==1
                alpha_k{1,id}=alpha_k{1,id}*(-1);
            else
                break;
            end
        end% line search spec==1
    end % end iter_line_search=1:ITER_MAX_LINE_SEARCH loop

    for i=1:n_var
        alpha_vec(1,i)=mean(alpha_k{i});
    end

end
