function [x_sol_cell,other_output_k,iter_info]=...
    Anderson_func(fun_fp,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.

n_var=size(x_0_cell,2);

spec=preliminary_spectral_func(spec,n_var);

ITER_MAX=spec.ITER_MAX;


t_Anderson=tic;

% varargin:1*XXX


FLAG_ERROR=0;

%% Read inputs
other_input_cell=varargin;


DIST_table=NaN(ITER_MAX,n_var);

for i=1:n_var
   elem_x(1,i)=prod(size(x_0_cell{i}));
end

resid_past_mat=NaN(sum(elem_x),ITER_MAX);
x_past_mat=NaN(sum(elem_x),ITER_MAX);
fun_k_past_mat=NaN(sum(elem_x),ITER_MAX);



%%%%%%%% Loop %%%%%%%%%%%
x_k_cell=x_0_cell;

FLAG_ERROR=0;
m=1;%%%%%

for k=0:ITER_MAX-1
    [fun_k_cell,other_output_k]=...
       fun_fp(x_k_cell{:},other_input_cell{:});

    loc=1;
    for i=1:n_var
        resid_past_mat(loc:loc+elem_x(1,i)-1,k+1)=fun_k_cell{i}(:)-x_k_cell{i}(:);
        x_past_mat(loc:loc+elem_x(1,i)-1,k+1)=x_k_cell{i}(:);
        fun_past_mat(loc:loc+elem_x(1,i)-1,k+1)=fun_k_cell{i}(:);
        loc=loc+elem_x(1,i);
    end

    if k>=1
        m_k=min(m,k);

        Y=resid_past_mat(:,k+1);%[]*1
        weight=1e-10;
        
        if 1==0
            %%% OLS spec 1
            X=Y-resid_past_mat(:,k-m_k+1:k);%[]*k
            X=X./max(max(abs(X)),1);
    
            alpha=(X'*X+weight*eye(size(X,2)))\(X'*Y);%k*1
    
            alpha_vec=[alpha;1-sum(alpha)];
        else

            %%% OLS spec 2
            X=diff(resid_past_mat(:,k-m_k+1:k+1),1,2);%[]*k
            X=X./max(max(abs(X)),1);
            weight=1e-0;
            gamma=(X'*X+weight*eye(size(X,2)))\(X'*Y);%k*1
    
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
        end % OLS spec


        beta_val=1;
        x_k_plus_1=beta_val*fun_past_mat(:,k-m_k+1:k+1)*alpha_vec+...
            (1-beta_val)*x_past_mat(k-m_k+1:k+1)*alpha_vec;%[]*1
        %x_k_plus_1=fun_past_mat(:,k);

        loc=1;
        for i=1:n_var
            x_k_plus_1_cell{1,i}=reshape(x_k_plus_1(loc:loc+elem_x(1,i)-1),size(x_0_cell{i}));
            loc=loc+elem_x(1,i);
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


end


