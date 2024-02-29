%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,n]=...
spectral_update_func(fun,x_k_cell,alpha_k,fun_k_cell,other_input_cell,...
n_var,line_search_spec,...
DIST_table,ITER_MAX_LINE_SEARCH,bound_spec,x_max_cell,x_min_cell)

    rho=0.8;
    M=10;gamma=10^(-4);

    %%% Update variables in the spectral algorithm %%%%%%%%%%%%%%%
    for n=1:ITER_MAX_LINE_SEARCH
     
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}-alpha_k{1,i}.*fun_k_cell{1,i};
       end % for loop wrt i

       if bound_spec==0
        [fun_k_plus_1_cell,other_output_k_plus_1]=...
           fun(x_k_plus_1_cell{:},other_input_cell{:});
       else
        [x_k_plus_1_cell,fun_k_plus_1_cell,other_output_k_plus_1]=...
            fun_bdd_func(...
            fun,x_k_plus_1_cell,other_input_cell,x_max_cell,x_min_cell,n_var);
       end

        %%% DIST: sup norm of F(x)=x-Phi(x). 
        DIST_vec=ones(1,n_var);
        for i=1:n_var
            DIST_vec(1,i)=max(abs(fun_k_plus_1_cell{1,i}),[],'all','omitnan');
        end
   
        if line_search_spec==1
            if n==3
                alpha_k{1,id}=alpha_k{1,id}*1;
            end

            DIST_PAST_MAX_vec=max(DIST_table(max(1,k+1-M+1):k+1,:),[],1);

            id=1;
            eta_k=DIST_table(1,id)/((1+k)^2);
            LHS=DIST_vec(id);

            %%%%% alpha_x^2+alpha_y^2???
            RHS=DIST_PAST_MAX_vec(id)+eta_k-gamma*(sum(alpha_k{1,id}(:).^2))*DIST_table(k+1,id);

            if LHS>RHS & mod(n,2)==0
                alpha_k{1,id}=alpha_k{1,id}.*(-rho);
            elseif LHS>RHS & mod(n,2)==1
                alpha_k{1,id}=alpha_k{1,id}*(-1);
            else
                break;
            end
        end% line search spec==1
    end % end n=1:ITER_MAX_LINE_SEARCH loop

end