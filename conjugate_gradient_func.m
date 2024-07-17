function d_k_cell=conjugate_gradient_func(d_k_cell,i,alpha_k,Delta_fun_cell,Delta_x_cell,fun_k_cell,fun_k_minus_1_cell)

                d_k_cell_original_i=d_k_cell{1,i};


               numer_i=reshape(alpha_k{1,i}.*Delta_fun_cell{1,i}-Delta_x_cell{1,i},[],1)'*...
                   fun_k_cell{1,i}(:);%%%%%
               numer_i=reshape(alpha_k{1,i}.*Delta_fun_cell{1,i},[],1)'*...
                   fun_k_cell{1,i}(:);
               %%numer_i=reshape(alpha_k{1,i}.*fun_k_cell{1,i},[],1)'*...
                   fun_k_cell{1,i}(:);

               %%% temp %%%
               %numer_i=alpha_k{1,i}*sqrt(sum(fun_k_cell{1,i}(:).^2));
               %denom_i=sqrt(sum(Delta_fun_cell{1,i}(:).^2));
               
               
               %denom_i=Delta_x_cell{1,i}(:)'*Delta_fun_cell{1,i}(:);%% Birgin Martinez (2001)
               denom_i=Delta_fun_cell{1,i}(:)'*Delta_fun_cell{1,i}(:);%%%% Faster ??              
               %%denom_i=fun_k_minus_1_cell{1,i}(:)'*fun_k_minus_1_cell{1,i}(:);             
               %%denom_i=Delta_fun_cell{1,i}(:)'*fun_k_minus_1_cell{1,i}(:);             
               
               %%% Dimension-wise??

               beta_k_i=numer_i./denom_i; %% Birgin and Martinez (2001)
               %%beta_k_i=min(max(beta_k_i,0),1);
               %%%beta_k_i
               
               d_k_cell{1,i}=d_k_cell{1,i}+beta_k_i.*Delta_x_cell{1,i};

               %%%%%%%%%%%%%%%%%%
               %if d_k_cell{1,i}(:)'*fun_k_cell{1,i}(:)>1e-3 %% cf. Birgin and Martinez (2001) eq (10)
               %    d_k_cell{1,i}=d_k_cell_original_i;
               %end
               %%%%%%%%%%%%%%%%%

end
