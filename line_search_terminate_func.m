function continue_backtracking_dummy=line_search_terminate_func(obj_val_vec,obj_val_table,d_k_cell,k,step_size,spec)

     M=spec.M;
     gamma=spec.gamma;
    
     LHS=sum(obj_val_vec);% ||F(x_{k+1})||^2

     obj_val_table_sum=sum(obj_val_table,2);
     obj_val_PAST_MAX=max(obj_val_table_sum(max(1,k+1-M+1):k+1,:));%max(||F(x_{k+1-M})||^2,...,||F(x_{k})||^2)
     eta_k=sqrt(obj_val_table_sum(1,:))/((1+k)^2);%obj_val: square of L2-norm
    
     f_k=obj_val_table_sum(k+1,:);%||F(x_{k})||^2
     RHS=obj_val_PAST_MAX+eta_k-(gamma*step_size^2)*(f_k);
   
    %%% If LHS<=RHS, exit the iteration; 
    %%% Otherwise (LHS>RHS), continue the iteration.
    continue_backtracking_dummy=(RHS-LHS<0);

    if 1==0
        %% Li Fukushima 
        d_k_L2_norm_squared=0;
        for i=1:size(d_k_cell,2)
            d_k_L2_norm_squared=d_k_L2_norm_squared+sum(d_k_cell{i}(:).^2);
        end

        LHS=sqrt(sum(obj_val_vec));%f(x_{k+1})
        RHS=(1+eta_k)*sqrt(f_k)-(gamma*step_size^2*d_k_L2_norm_squared);

        continue_backtracking_dummy=(RHS-LHS<0);
    end
end
