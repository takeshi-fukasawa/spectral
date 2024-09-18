function continue_backtracking_dummy=line_search_terminate_func(obj_val_vec,obj_val_table,n_var,k,step_size,spec)

     M=spec.M;
     gamma=spec.gamma;
     
     obj_val_table_sum=sum(obj_val_table,2);
     LHS=sum(obj_val_vec);%1*n_var
       
     obj_val_PAST_MAX=max(obj_val_table_sum(max(1,k+1-M+1):k+1,:));
     eta_k=sqrt(obj_val_table_sum(1,:))/((1+k)^2);
    
     eta_k=1*(obj_val_table_sum(1,:))/((1+k)^2);%%%%
    
     RHS=obj_val_PAST_MAX+eta_k-(gamma*step_size^2)*(obj_val_table_sum(k+1,:));

         
    %%% If LHS<=RHS, exit the iteration; 
    %%% Otherwise (LHS>RHS), continue the iteration.
    continue_backtracking_dummy=(RHS-LHS<0);


end
