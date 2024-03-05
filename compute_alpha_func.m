function alpha_k_i=compute_alpha_func(sum_Delta_x_fun_i,...
    sum_Delta_x_x_i,sum_Delta_fun_fun_i,...
    compute_alpha_spec,dampening_param)

if compute_alpha_spec>=3
    %%%% Specification proposed in Varadhan and Roland (2008)
    alpha_k_i=sqrt(sum_Delta_x_x_i./sum_Delta_fun_fun_i);%scalar or vector (wrt the dimension specified in "update_spec")

    if compute_alpha_spec==4
        %%% Sign can be negative
        alpha_k_i=alpha_k_i.*sign(sum_Delta_x_fun_i);
    end

elseif compute_alpha_spec==1
    %%%% BB first spec
    alpha_k_i=sum_Delta_x_fun_i./sum_Delta_fun_fun_i;%scalar or vector (wrt the dimension specified in "update_spec")

elseif compute_alpha_spec==2
    %%%% BB second spec
    alpha_k_i=sum_Delta_x_x_i./sum_Delta_x_fun_i;%scalar or vector (wrt the dimension specified in "update_spec")
end % compute_alpha_spec

    alpha_k_i((isnan(alpha_k_i)==1))=1;%%%
    alpha_k_i((isinf(alpha_k_i)==1))=1;%%%
    alpha_k_i(((alpha_k_i==0)))=1;%%%

    
    if isempty(dampening_param)==0
        alpha_k_i=alpha_k_i*dampening_param;
    end

end
