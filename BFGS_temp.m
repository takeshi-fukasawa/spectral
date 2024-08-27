%%%%%%%% Loop %%%%%%%%%%%

if conv==0 & ITER_MAX>=2
for k=0:ITER_MAX-2


   if k>=1

      for i=1:n_var
        Delta_x_cell{1,i}=x_k_cell{i}-x_k_minus_1_cell{i};
        Delta_fun_cell{1,i}=fun_k_cell{i}-fun_k_minus_1_cell{i};
      end % loop wrt i

      if spec.update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;
          end

      elseif spec.BFGS_spec==0
        [alpha_k,alpha_max]=compute_alpha_func(...
         Delta_x_cell,Delta_fun_cell,spec,k,DIST_table(k+1,:));
         spec.alpha_max=alpha_max;

      elseif spec.BFGS_spec==1
         s=Delta_x_cell{1,1}(:);

         if sum(abs(s(:)))==0 & spec.dampening_param{1,1}>0
            H_k=H_k_minus_1;
         else
            y=Delta_fun_cell{1,1}(:);
            
            if spec.fixed_point_iter_spec==1
                y=y.*(-1);%%%####
            end

            I=eye(size(s,1));
            syp=s*y';
            spy=s'*y;% scalar
            ysp=y*s';
            ssp=s*s';

            H_k=(I-syp./spy)*H_k_minus_1*(I-ysp./spy)+ssp./spy;
        end

    
    end

  else % k==0

      if spec.BFGS_spec==1
         H_k=eye(size(fun_k_cell{1,1}(:),1));
      end
   end

    %%% Update variables %%%%%%%%%%%%%%%
   for i=1:n_var
        

        if spec.BFGS_spec==1 & k>=2 & i==1

            if alpha_k{1}==0
                H_k=eye(size(d_k_cell{1,i}(:),1));
            end

            d_k_cell{1,1}=H_k*d_k_cell{1,1}(:);
        end

        alpha_table(k+1,i)=max(alpha_k{1,i}(:));
    end % for loop wrt i

    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,obj_val_vec,iter_line_search,step_size]=...
        spectral_update_func(fun,x_k_cell,fun_k_cell,d_k_cell,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);

        

    %%% Replace variables for the next iteration
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;
    fun_k_cell=fun_k_plus_1_cell;
   

    if spec.BFGS_spec==1
       H_k_minus_1=H_k;

       if isnan(sum(H_k_minus_1(:)))==1
            temp=0;
            temp2=0;
       end
    end


end %% end of for loop wrt k=0:ITER_MAX-1

