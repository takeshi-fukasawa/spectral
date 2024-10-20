function [x_sol_cell,other_output_k,iter_info]=...
    Anderson_func(fun_fp,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.

warning('off','all')

n_var=size(x_0_cell,2);

spec=preliminary_spectral_func(spec,n_var);

ITER_MAX=spec.ITER_MAX;

if isfield(spec,'m_Anderson')==0
    spec.m=5;
else
    spec.m=spec.m_Anderson;
end

if isfield(spec,'type_Anderson')==0
    spec.type_Anderson=2;% Type II Anderson
end

m=spec.m;
type_Anderson=spec.type_Anderson;



% varargin:1*XXX


FLAG_ERROR=0;

%% Read inputs
other_input_cell=varargin;


DIST_table=NaN(ITER_MAX,1);

for i=1:n_var
   elem_x(1,i)=prod(size(x_0_cell{i}));
end

t_Anderson=tic;



%%%%%%%% Loop %%%%%%%%%%%

FLAG_ERROR=0;

x_k_cell=x_0_cell;

for k=0:ITER_MAX-1

    [fun_k_cell,other_output_k]=...
       fun_fp(x_k_cell{:},other_input_cell{:});

   fun_k_vec=cell_to_vec_func(fun_k_cell);
   x_k_vec=cell_to_vec_func(x_k_cell);
   resid_k_vec=fun_k_vec-x_k_vec;

   if k==0
      resid_past_mat=resid_k_vec;
      fun_past_mat=fun_k_vec;
      x_past_mat=x_k_vec;
   else
     resid_past_mat=[resid_past_mat,resid_k_vec];%[]*(k+1)
     x_past_mat=[x_past_mat,x_k_vec];
     fun_past_mat=[fun_past_mat,fun_k_vec];
   end 

    if k>=1
        m_k=min(m,k);

        Z=resid_past_mat(:,k+1);%[]*1
        


        DF=diff(resid_past_mat(:,k-m_k+1:k+1),1,2);%[]*k
        %%DF=DF./max(max(abs(DF)),1);

            if isnan(sum(DF(:)))==1
                FLAG_ERROR=1;
                break;
            end

            if type_Anderson==1
                %%% Type I Anderson (Corresponding to Good Broyden update)%%% 
                DX=diff(x_past_mat(:,k-m_k+1:k+1),1,2);%[]*k
                %DX=DX./max(max(abs(DX(:))),1);%%%
                gamma=(DX'*DF)\(DX'*Z);
            
            else % type_Anderson==2
                %%% Type II Anderson (Corresponding to Good Broyden update)
                %%% Based on the computation of DF'*DF
                gamma=(DF'*DF)\(DF'*Z);%k*1
                

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
            end


        beta_val=1;
       if beta_val==1
            x_k_plus_1=fun_past_mat(:,k-m_k+1:k+1)*alpha_vec;%[]*1
       else
            x_k_plus_1=beta_val*fun_past_mat(:,k-m_k+1:k+1)*alpha_vec+...
                (1-beta_val)*x_past_mat(k-m_k+1:k+1)*alpha_vec;%[]*1
       end
       %x_k_plus_1=fun_past_mat(:,k);

       x_k_plus_1_cell=vec_to_cell_func(x_k_plus_1,elem_x,x_0_cell);

        if 1==0 % Stabilization ??
            [fun_k_plus_1_cell,other_output_k]=...
            fun_fp(x_k_plus_1_cell{:},other_input_cell{:});
    
            DIST_fp=0;DIST_Anderson=0;
            for i=1:n_var
                DIST_fp=DIST_fp+sum(fun_k_cell{i}(:)-x_k_cell{i}(:));
                DIST_Anderson=DIST_Anderson+sum(fun_k_plus_1_cell{i}(:)-x_k_plus_1_cell{i}(:));
            end
    
            if DIST_fp<DIST_Anderson
                x_k_plus_1_cell=fun_k_cell;
            end
        end


           
        DIST=max(abs(resid_past_mat(:,k+1)));

        interval=10;
        if k-floor(k/interval)*interval==0&spec.DEBUG==1
            DIST
        end

        if DIST<spec.TOL
            %DIST
            FLAG_ERROR=0;
            break;
        end
        DIST_table(k)=DIST;
        
    else % k==0
        x_k_plus_1_cell=fun_k_cell;
    end % k==0 or not

    x_k_cell=x_k_plus_1_cell;



end %% end of for loop wrt k=0:ITER_MAX-1


%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc(t_Anderson);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=iter_info.n_iter;

iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

iter_info.DIST_table=DIST_table;
iter_info.type_Anderson=type_Anderson;


end


