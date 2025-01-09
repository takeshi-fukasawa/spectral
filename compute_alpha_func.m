function [alpha_k,alpha_max]=compute_alpha_func(...
    Delta_x_cell,Delta_fun_cell,spec,k)
    
    n_var=size(Delta_x_cell,2);
    
    for i=1:n_var
          if isempty(spec.dim_hetero_alpha)==1 || spec.dim_hetero_alpha(i)==0
            sum_dim_ids='all';
    
          else  %%% XXX-dependent tune parameters
            sum_dim_ids=1:size(size(Delta_x_cell{1,i}),2);
            sum_dim_ids=sum_dim_ids(sum_dim_ids~=(spec.dim_hetero_alpha(i)));
          end 
    
            sum_Delta_x_x_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_x_cell{1,i},sum_dim_ids,'omitnan');%vector
            sum_Delta_fun_fun_cell{1,i}=sum(Delta_fun_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector      
            sum_Delta_x_fun_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector  
    end % loop wrt i
    
    if spec.common_spectral_coef_spec==1
        sum_Delta_x_x=0;
        sum_Delta_fun_fun=0;
        sum_Delta_x_fun=0;

        i_start=1;

        for i=i_start:n_var
            sum_Delta_x_x=sum_Delta_x_x+sum_Delta_x_x_cell{i};
            sum_Delta_fun_fun=sum_Delta_fun_fun+sum_Delta_fun_fun_cell{i};
            sum_Delta_x_fun=sum_Delta_x_fun+sum_Delta_x_fun_cell{i};
        end
    
        for i=i_start:n_var
            sum_Delta_x_x_cell{1,i}=sum_Delta_x_x;
            sum_Delta_fun_fun_cell{1,i}=sum_Delta_fun_fun;
            sum_Delta_x_fun_cell{1,i}=sum_Delta_x_fun;
        end % for loop wrt i
    end % if statement
    
    for i=1:n_var
    if spec.spectral_coef_spec>=3
        %%%% Specification proposed in Varadhan and Roland (2008)
        alpha_k{1,i}=sqrt(sum_Delta_x_x_cell{1,i}./sum_Delta_fun_fun_cell{1,i});%scalar or vector (wrt the dimension specified in "dim_hetero_alpha")
    
        if spec.spectral_coef_spec==4
            %%% Sign can be negative
            alpha_k{1,i}=-alpha_k{1,i}.*sign(sum_Delta_x_fun_cell{1,i});
        end % if statement
    
    elseif spec.spectral_coef_spec==1
        %%%% Barzilai and Borwein (1988) first spec
        alpha_k{1,i}=-sum_Delta_x_fun_cell{1,i}./sum_Delta_fun_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "dim_hetero_alpha")
    
    elseif spec.spectral_coef_spec==2
        %%%% Barzilai and Borwein (1988) second spec
        alpha_k{1,i}=-sum_Delta_x_x_cell{1,i}./sum_Delta_x_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "dim_hetero_alpha")
    end % spectral_coef_spec
        
        alpha_k{1,i}((isnan(alpha_k{1,i})==1))=1;%%%
        alpha_k{1,i}((isinf(alpha_k{1,i})==1))=1;%%%
        alpha_k{1,i}(((alpha_k{1,i}==0)))=1;%%%
    
        alpha_max=spec.alpha_max;
        alpha_k{1,i}=min(alpha_max,alpha_k{1,i});%%%%%
    
        alpha_min=spec.alpha_min;
        alpha_k{1,i}=max(alpha_min,alpha_k{1,i});%%%%%
    
        if isempty(spec.dampening_param)==0
           dampening_param=spec.dampening_param;
            alpha_k{1,i}=alpha_k{1,i}*dampening_param{1,i};
        end
    
       end % for loop wrt i
    
    end
    