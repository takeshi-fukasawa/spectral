function continue_backtracking_dummy=line_search_terminate_func(obj_val_vec,obj_val_table,n_var,k,step_size,spec)

     M=spec.M;
     gamma=spec.gamma;
     
     LHS=obj_val_vec;%1*n_var
           
     for i=1:n_var 
         obj_val_PAST_MAX_i=max(obj_val_table(max(1,k+1-M+1):k+1,i));

         eta_k_i=sqrt(obj_val_table(1,i))/((1+k)^2);
         RHS_i=obj_val_PAST_MAX_i+eta_k_i-(gamma*step_size^2)*(obj_val_table(k+1,i));%1*n_var
              
         
         %%% If LHS<=RHS, exit the iteration; 
         %%% Otherwise (LHS>RHS), continue the iteration.
         continue_backtracking_dummy(1,i)=(RHS_i-LHS(1,i)<0);

         %if sum(continue_backtracking_dummy(:)==n_var)
         %   k
         %   RHS_i-LHS(1,i)
         %   continue_backtracking_dummy(1,i)
         %end


     end % for loop wrt i




end
