function [alpha_k,alpha_max]=compute_alpha_func(...
Delta_x_cell,Delta_fun_cell,spec,k,norm_fun_k)

n_var=size(Delta_x_cell,2);

update_spec=spec.update_spec;

for i=1:n_var
      if isempty(update_spec)==1
        sum_dim_ids='all';

      elseif sum(update_spec(:))>0 %%% XXX-dependent tune parameters
         sum_dim_ids=1:size(size(Delta_x_cell{1,i}),2);
                  sum_dim_ids=sum_dim_ids(sum_dim_ids~=update_spec);
         
      end % isempty(update_spec)==1 or others?

      %%if isempty(update_spec)==1 | (isempty(update_spec)==0 & sum(update_spec(:))>0) 
        sum_Delta_x_x_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_x_cell{1,i},sum_dim_ids,'omitnan');%vector
        sum_Delta_fun_fun_cell{1,i}=sum(Delta_fun_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector      
        sum_Delta_x_fun_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector  
      %%end % if statement

end % loop wrt i

if spec.common_alpha_spec==1
    sum_Delta_x_x=0;
    sum_Delta_fun_fun=0;
    sum_Delta_x_fun=0;
    for i=1:n_var
        sum_Delta_x_x=sum_Delta_x_x+sum_Delta_x_x_cell{i};
        sum_Delta_fun_fun=sum_Delta_fun_fun+sum_Delta_fun_fun_cell{i};
        sum_Delta_x_fun=sum_Delta_x_fun+sum_Delta_x_fun_cell{i};
    end

for i=1:n_var
    sum_Delta_x_x_cell{1,i}=sum_Delta_x_x;
    sum_Delta_fun_fun_cell{1,i}=sum_Delta_fun_fun;
    sum_Delta_x_fun_cell{1,i}=sum_Delta_x_fun;
end % for loop wrt i
     end % if statement

for i=1:n_var
if spec.compute_alpha_spec>=3
    %%%% Specification proposed in Varadhan and Roland (2008)
    alpha_k{1,i}=sqrt(sum_Delta_x_x_cell{1,i}./sum_Delta_fun_fun_cell{1,i});%scalar or vector (wrt the dimension specified in "update_spec")

    if spec.compute_alpha_spec==4
        %%% Sign can be negative
        alpha_k{1,i}=-alpha_k{1,i}.*sign(sum_Delta_x_fun_cell{1,i});
    end % if statement

elseif spec.compute_alpha_spec==1
    %%%% BB first spec
    alpha_k{1,i}=-sum_Delta_x_fun_cell{1,i}./sum_Delta_fun_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "update_spec")

elseif spec.compute_alpha_spec==2
    %%%% BB second spec
    alpha_k{1,i}=-sum_Delta_x_x_cell{1,i}./sum_Delta_x_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "update_spec")
end % compute_alpha_spec
    
    alpha_k{1,i}((isnan(alpha_k{1,i})==1))=1;%%%
    alpha_k{1,i}((isinf(alpha_k{1,i})==1))=1;%%%
    alpha_k{1,i}(((alpha_k{1,i}==0)))=1;%%%

    %%%alpha_max=10;%%%%
    alpha_max=spec.alpha_max;
    alpha_k{1,i}=min(alpha_max,alpha_k{1,i});%%%%%

    if spec.stabilize_spec==1 & 1==0
      Delta=spec.Delta;
      Delta=0.01;%%%
	    alpha_k{1,i}=min(alpha_k{1,i},Delta./norm_fun_k(1,i));
    end

    %%%%%%%
    %%%if k<=300
    %%%    alpha_k{1,i}=1;%%%%
    %%%end

    %temp=(alpha_k{1,i}==alpha_max);
    %if sum(temp(:))>0
    %    if alpha_max<=5
    %        step_factor=1.5;
    %        alpha_max=alpha_max*step_factor;
    %    end
    %end
    %%%%%%

    if isempty(spec.dampening_param)==0
       dampening_param=spec.dampening_param;
        alpha_k{1,i}=alpha_k{1,i}*dampening_param{1,i};
    end

   end % for loop wrt i

end
