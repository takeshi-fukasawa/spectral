function continue_backtracking_dummy=line_search_terminate_func(obj_val_vec,obj_val_table,n_var,k,step_size,fun_k_cell,d_k_cell,spec)

     M=spec.M;
     gamma=spec.gamma;
     
     LHS=obj_val_vec;%1*n_var
           
     for i=1:n_var 
         obj_val_PAST_MAX_i=max(obj_val_table(max(1,k+1-M+1):k+1,i));

         if spec.nonlinear_eq_spec(1,i)==1   %% Nonlinear eq; La Cruz (2006)      
              eta_k_i=sqrt(obj_val_table(1,i))/((1+k)^2);
              RHS_i=obj_val_PAST_MAX_i+eta_k_i-(gamma*step_size^2)*(obj_val_table(k+1,i));%1*n_var
              
         elseif spec.BFGS_spec==0 % minimization problem
              RHS_i=obj_val_PAST_MAX_i+gamma*step_size*fun_k_cell{1,i}(:)'*d_k_cell{1,i}(:);%1*n_var

         elseif spec.BFGS_spec==1 % minimization problem
              RHS_i=obj_val_table(k+1,i)+gamma*step_size*fun_k_cell{1,i}(:)'*d_k_cell{1,i}(:);%1*n_var

              %%eta_k_i=sqrt(obj_val_table(1,i))/((1+k)^2);
              %%RHS_i=obj_val_PAST_MAX_i+eta_k_i-(gamma*step_size^2)*(obj_val_table(k+1,i));%1*n_var%%%%%%%=> not converge....? (BLP)
              %%%RHS_i=obj_val_PAST_MAX_i+gamma*step_size*fun_k_cell{1,i}(:)'*d_k_cell{1,i}(:);%1*n_var


         end %if spec.minimization_spec(1,i)==0

         %%% If LHS<=RHS, exit the iteration; 
         %%% Otherwise (LHS>RHS), continue the iteration.
         continue_backtracking_dummy(1,i)=(RHS_i-LHS(1,i)<0);

     end % for loop wrt i

     if 1==1
          weight_param=10;
          %%% Structural est
          RHS=obj_val_table(k+1,1)+weight_param*sum(obj_val_table(k+1,2:end))+gamma*step_size*fun_k_cell{1,1}(:)'*d_k_cell{1,1}(:);%1*n_var
          LHS=obj_val_vec(1,1)+weight_param*sum(obj_val_vec(1,2:end));
          continue_backtracking_dummy(1,1:n_var)=(RHS-LHS<0);

          %%[LHS,RHS]
     end


end
