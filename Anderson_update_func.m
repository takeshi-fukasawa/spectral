function [x_k_plus_1,FLAG_ERROR]=Anderson_update_func(resid_past_mat,x_past_mat,fun_fp_past_mat,type_Anderson,m,k,beta_val,FLAG_ERROR)

        m_k=min(m,k);
 
           Z=resid_past_mat(:,k+1);%[]*1; f_{k}
        
           DF=diff(resid_past_mat(:,k-m_k+1:k+1),1,2);%[]*m_k;Df_{k-mk},...,Df_{k-1}; Df_{k-1}=f_{k}-f_{k-1}

            if isnan(sum(DF(:)))==1
                warning ("NaN DF")
                FLAG_ERROR=1;
            end

            if type_Anderson==1
                %%% Type I Anderson (Corresponding to Good Broyden update)%%% 
                DX=diff(x_past_mat(:,k-m_k+1:k+1),1,2);%[]*m_k; Dx_{k-mk},...,Dx_{k-1}
                gamma=(DX'*DF)\(DX'*Z);
            
            else % type_Anderson==2
                %%% Type II Anderson (Corresponding to Good Broyden update)
                %%% Based on the computation of DF'*DF
                gamma=(DF'*DF)\(DF'*Z);%m_k*1
                

                %%% Based on QR factorization (Slow??)
                %[Q,R]=qr(DF);
                %Q1=Q(:,1:size(DF,2));
                %R1=R(1:size(DF,2),:);
                %gamma=inv(R1)*Q1'*Z;

                %%% Based on SVD
                %[U,S,V] = svd(DF);
                %U1=U(:,1:size(X,2));
                %S1=S(1:size(X,2),:);
                %gamma=V*inv(S1)*(U1')*Z;
            end%Type_I_Anderson_spec==1 or 0??


            alpha_vec=zeros(size(gamma,1)+1,1);
            for id=1:size(alpha_vec,1)
                if id==1
                    alpha_vec(1)=gamma(1);
                elseif id>=2 & id<=size(alpha_vec,1)-1
                    alpha_vec(id)=gamma(id)-gamma(id-1);
                elseif id==size(alpha_vec,1)
                    alpha_vec(id)=1-gamma(id-1);
                end
            end%id


     %% Update x_k_plus_1

       if beta_val==1
            x_k_plus_1=fun_fp_past_mat(:,k-m_k+1:k+1)*alpha_vec;%[]*1
       else
            x_k_plus_1_i=beta_val*fun_fp_past_mat(:,k-m_k+1:k+1)*alpha_vec+...
                (1-beta_val)*x_past_mat(k-m_k+1:k+1)*alpha_vec;%[]*1
       end%beta_val==1
end

