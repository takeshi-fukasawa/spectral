function [alpha_k,alpha_max]=compute_alpha_func(...
Delta_x_cell,Delta_fun_cell,...
    common_alpha_spec,compute_alpha_spec,dampening_param,update_spec,...
    alpha_max,k)

n_var=size(Delta_x_cell,2);


for i=1:n_var
      if isempty(update_spec)==1
        sum_dim_ids='all';

      elseif sum(update_spec(:))>0 %%% XXX-dependent tune parameters
         sum_dim_ids=1:size(size(Delta_x_cell{1,i}),2);
         %%%sum_dim_ids=sum_dim_ids(sum_dim_ids~=update_spec(i));
         sum_dim_ids=sum_dim_ids(sum_dim_ids~=update_spec);
         
      end % isempty(update_spec)==1 or others?

      %%if isempty(update_spec)==1 | (isempty(update_spec)==0 & sum(update_spec(:))>0) 
        sum_Delta_x_x_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_x_cell{1,i},sum_dim_ids,'omitnan');%vector
        sum_Delta_fun_fun_cell{1,i}=sum(Delta_fun_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector      
        sum_Delta_x_fun_cell{1,i}=sum(Delta_x_cell{1,i}.*Delta_fun_cell{1,i},sum_dim_ids,'omitnan');%vector  
      %%end % if statement

end % loop wrt i

if common_alpha_spec==1
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
if compute_alpha_spec>=3
    %%%% Specification proposed in Varadhan and Roland (2008)
    alpha_k{1,i}=sqrt(sum_Delta_x_x_cell{1,i}./sum_Delta_fun_fun_cell{1,i});%scalar or vector (wrt the dimension specified in "update_spec")

    if compute_alpha_spec==4
        %%% Sign can be negative
        alpha_k{1,i}=-alpha_k{1,i}.*sign(sum_Delta_x_fun_cell{1,i});
    end % if statement

elseif compute_alpha_spec==1
    %%%% BB first spec
    alpha_k{1,i}=-sum_Delta_x_fun_cell{1,i}./sum_Delta_fun_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "update_spec")

elseif compute_alpha_spec==2
    %%%% BB second spec
    alpha_k{1,i}=-sum_Delta_x_x_cell{1,i}./sum_Delta_x_fun_cell{1,i};%scalar or vector (wrt the dimension specified in "update_spec")
end % compute_alpha_spec
    
    alpha_k{1,i}((isnan(alpha_k{1,i})==1))=1;%%%
    alpha_k{1,i}((isinf(alpha_k{1,i})==1))=1;%%%
    alpha_k{1,i}(((alpha_k{1,i}==0)))=1;%%%

    %%%alpha_max=10;%%%%
    alpha_k{1,i}=min(alpha_max,alpha_k{1,i});%%%%%

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

    if isempty(dampening_param)==0
        alpha_k{1,i}=alpha_k{1,i}*dampening_param{1,i};
    end

   end % for loop wrt i

end
