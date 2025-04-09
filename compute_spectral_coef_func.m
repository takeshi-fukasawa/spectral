function [spectral_coef_k,spectral_coef_max]=compute_spectral_coef_func(...
    Delta_x_cell,Delta_fun_cell,spec,k)
    
    n_var=size(Delta_x_cell,2);
    
    for i=1:n_var
          if isempty(spec.dim_hetero_spectral_coef)==1 || spec.dim_hetero_spectral_coef(i)==0
            sum_dim_ids='all';
    
          else  %%% XXX-dependent tune parameters
            sum_dim_ids=1:size(size(Delta_x_cell{1,i}),2);
            sum_dim_ids=sum_dim_ids(sum_dim_ids~=(spec.dim_hetero_spectral_coef(i)));
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
        spectral_coef_k{1,i}=sqrt(sum_Delta_x_x_cell{1,i}./sum_Delta_fun_fun_cell{1,i});%scalar or vector (wrt the dimension specified in "dim_hetero_spectral_coef")
    
        if spec.spectral_coef_spec==4
            %%% Sign can be negative
            spectral_coef_k{1,i}=-spectral_coef_k{1,i}.*sign(sum_Delta_x_fun_cell{1,i});
        end % if statement
    
    elseif spec.spectral_coef_spec==1
        %%%% Barzilai and Borwein (1988) first spec
        spectral_coef_k{1,i}=-sum_Delta_x_fun_cell{1,i}./sum_Delta_fun_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "dim_hetero_spectral_coef")
    
    elseif spec.spectral_coef_spec==2
        %%%% Barzilai and Borwein (1988) second spec
        spectral_coef_k{1,i}=-sum_Delta_x_x_cell{1,i}./sum_Delta_x_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "dim_hetero_spectral_coef")
    end % spectral_coef_spec
        
        spectral_coef_k{1,i}((isnan(spectral_coef_k{1,i})==1))=1;%%%
        spectral_coef_k{1,i}((isinf(spectral_coef_k{1,i})==1))=1;%%%
        spectral_coef_k{1,i}(((spectral_coef_k{1,i}==0)))=1;%%%
    
        spectral_coef_max=spec.spectral_coef_max;
        spectral_coef_k{1,i}=min(spectral_coef_max,spectral_coef_k{1,i});%%%%%
    
        spectral_coef_min=spec.spectral_coef_min;
        spectral_coef_k{1,i}=max(spectral_coef_min,spectral_coef_k{1,i});%%%%%
    
        if isempty(spec.dampening_param)==0
           dampening_param=spec.dampening_param;
            spectral_coef_k{1,i}=spectral_coef_k{1,i}*dampening_param{1,i};
        end
    
       end % for loop wrt i
    
    end
    