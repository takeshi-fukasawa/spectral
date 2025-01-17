function [x_sol_cell,other_output_k,iter_info]=...
    Anderson_func2(fun_fp,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in January 2025.
%%% Variable-type-specific updates introduced

warning('off','all')

n_var=size(x_0_cell,2);

spec=preliminary_setting_func(spec,n_var);

x_max_cell=spec.x_max_cell;
x_min_cell=spec.x_min_cell;

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

for i=1:n_var
   elem_x(1,i)=prod(size(x_0_cell{i}));
end

t_Anderson=tic;

if spec.common_Anderson_coef_spec==1
  n_var_type=1;
else
   n_var_type=n_var;
end
DIST_table=NaN(ITER_MAX,n_var_type);
DIST_vec=NaN(1,n_var_type);


%%%%%%%% Loop %%%%%%%%%%%

FLAG_ERROR=0;

x_k_cell=x_0_cell;
x_k_plus_1_cell=[];
d_k_cell=[];

[fun_fp_k_cell,other_output_k]=...
       fun_fp(x_k_cell{:},other_input_cell{:});

fun_fp_k_plus_1_cell=fun_fp_k_cell;

for i=1:n_var_type
     DIST_vec(i)=norm_func(fun_fp_k_cell{i}(:)-x_k_cell{i}(:),x_k_cell{i}(:),spec.norm_spec(i));
end

obj_val_table=NaN(ITER_MAX,n_var_type);
obj_val_table(1,:)=DIST_vec.^2;% L2 norm

step_size_table=NaN(ITER_MAX,n_var_type);

for k=0:ITER_MAX-1

    if spec.common_Anderson_coef_spec==1
        x_k_vec_cell{1}=cell_to_vec_func(x_k_cell);
        fun_fp_k_vec_cell{1}=cell_to_vec_func(fun_fp_k_cell);
    else %spec.common_Anderson_coef_spec==[]
        for i=1:n_var
            x_k_vec_cell{i}=x_k_cell{i}(:);
            fun_fp_k_vec_cell{i}=fun_fp_k_cell{i}(:);
       end
   end

    for i=1:n_var_type
       resid_k_vec_cell{i}=fun_fp_k_vec_cell{i}-x_k_vec_cell{i};
    end
    

   if k==0
        for i=1:n_var_type
          resid_past_mat_cell{i}=resid_k_vec_cell{i};
          fun_fp_past_mat_cell{i}=fun_fp_k_vec_cell{i};
          x_past_mat_cell{i}=x_k_vec_cell{i};
       end

   else%k>=1

        for i=1:n_var_type
          resid_past_mat_cell{i}=[resid_past_mat_cell{i},resid_k_vec_cell{i}];%[]*(k+1)
          fun_fp_past_mat_cell{i}=[fun_fp_past_mat_cell{i},fun_fp_k_vec_cell{i}];
          x_past_mat_cell{i}=[x_past_mat_cell{i},x_k_vec_cell{i}];

          
        end
   end 

    if k>=1
        m_k=min(m,k);

        for i=1:n_var_type 
            Z=resid_past_mat_cell{i}(:,k+1);%[]*1
        
           DF=diff(resid_past_mat_cell{i}(:,k-m_k+1:k+1),1,2);%[]*k
           %%DF=DF./max(max(abs(DF)),1);

            if isnan(sum(DF(:)))==1
                FLAG_ERROR=1;
                break;
            end

            if type_Anderson==1
                %%% Type I Anderson (Corresponding to Good Broyden update)%%% 
                DX=diff(x_past_mat_cell{i}(:,k-m_k+1:k+1),1,2);%[]*k
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
            end%id


     %% Update x_k_plus_1
        beta_val=1;
       if beta_val==1
            x_k_plus_1_i=fun_fp_past_mat_cell{i}(:,k-m_k+1:k+1)*alpha_vec;%[]*1
       else
            x_k_plus_1_i=beta_val*fun_fp_past_mat_cell{i}(:,k-m_k+1:k+1)*alpha_vec+...
                (1-beta_val)*x_past_mat_cell{i}(k-m_k+1:k+1)*alpha_vec;%[]*1
       end%beta_val==1

       if spec.common_Anderson_coef_spec==1
            x_k_plus_1_cell=vec_to_cell_func(x_k_plus_1_i,elem_x,x_0_cell);
            d_k_cell={x_k_plus_1_cell{1}-x_k_cell{1}};
       else
           x_k_plus_1_cell{i}=reshape(x_k_plus_1_i,size(x_0_cell{i}));
           d_k_cell{i}={x_k_plus_1_cell{i}-x_k_cell{i}};
       end

    end%i=1:n_var_type

    if spec.line_search_spec==0
       [fun_fp_k_plus_1_cell,other_output_k_plus_1]=...
       fun_fp(x_k_plus_1_cell{:},other_input_cell{:});

       for i=1:n_var_type
          DIST_vec(i)=norm_func(fun_fp_k_plus_1_cell{i}(:)-x_k_plus_1_cell{i}(:),reshape(x_past_mat_cell{i}(:,1),[],1),spec.norm_spec(i));
       end
      step_size=ones(1,n_var_type);
      obj_val_vec=DIST_vec.^2;


    else%spec.line_search_spec==1
        [x_k_plus_1_cell, fun_k_plus_1_cell,...
        other_output_k_plus_1,DIST_vec,iter_line_search,...
        obj_val_vec,step_size]=...
        update_func(fun_fp,x_k_cell,d_k_cell,other_input_cell,...
        n_var_type,spec,x_max_cell,x_min_cell,k,obj_val_table);

        for i=1:n_var_type
            fun_fp_k_plus_1_cell{i}=fun_k_plus_1_cell{i}+x_k_plus_1_cell{i};
        end
    end%spec.line_search_spec==0 or 1


    other_output_k=other_output_k_plus_1;

        interval=10;
        if k-floor(k/interval)*interval==0&spec.DEBUG==1
            DIST_vec
        end
        DIST=max(DIST_vec);

        conv=(sum((DIST_vec<spec.TOL),'all')==n_var_type);
        if conv==1
            FLAG_ERROR=0;
            break;
        end
        DIST_table(k,:)=DIST_vec;
        obj_val_table(k,:)=obj_val_vec;
        step_size_table(k,:)=step_size;
        
    else % k==0
        x_k_plus_1_cell=fun_fp_k_cell;
    end % k==0 or not

    x_k_cell=x_k_plus_1_cell;
    fun_fp_k_cell=fun_fp_k_plus_1_cell;

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
iter_info.step_size_table=step_size_table;
iter_info.type_Anderson=type_Anderson;
iter_info.spec=spec;

end


