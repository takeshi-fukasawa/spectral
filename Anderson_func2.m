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

if spec.common_Anderson_coef_spec==1
  n_var_type=1;
else
   n_var_type=n_var;
end

DIST_table=NaN(ITER_MAX,n_var);
DIST_vec=NaN(1,n_var_type);
obj_val_vec=NaN(1,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);
obj_val_table=NaN(ITER_MAX,n_var);
step_size_table=NaN(ITER_MAX,n_var);

FLAG_ERROR=0;
%%%%%%%% Iteration %%%%%%%%%%%
t_Anderson=tic;


x_k_cell=x_0_cell;
x_k_plus_1_cell=[];
d_k_cell=[];

feval=0;

%%% k==0
k=0;
[fun_fp_k_cell,other_output_k]=...
   fun_fp(x_k_cell{:},other_input_cell{:});

feval=feval+1;

for i=1:n_var
     DIST_vec(i)=norm_func(fun_fp_k_cell{i}(:)-x_k_cell{i}(:),x_k_cell{i}(:),spec.norm_spec(i));
     obj_val_table(1,i)=sum((fun_fp_k_cell{i}(:)-x_k_cell{i}(:)).^2);%(L2 norm)^2
end

DIST_table(1,:)=DIST_vec;

if max(DIST_vec)>=spec.TOL & ITER_MAX>=2
 
    if spec.common_Anderson_coef_spec==1
        x_k_vec_cell{1}=cell_to_vec_func(x_k_cell);
        fun_fp_k_vec_cell{1}=cell_to_vec_func(fun_fp_k_cell);
    else %spec.common_Anderson_coef_spec==[]
        for i=1:n_var
            x_k_vec_cell{i}=x_k_cell{i}(:);
            fun_fp_k_vec_cell{i}=fun_fp_k_cell{i}(:);
       end
   end% if spec.common_Anderson_coef_spec==1 or 0

        for i=1:n_var_type
          fun_fp_past_mat_cell{i}=fun_fp_k_vec_cell{i};
          x_past_mat_cell{i}=x_k_vec_cell{i};
          resid_past_mat_cell{i}=fun_fp_past_mat_cell{i}-x_past_mat_cell{i};
       end

%%% k==1
k=1;
x_k_cell=fun_fp_k_cell;
   [fun_fp_k_cell,other_output_k]=...
       fun_fp(x_k_cell{:},other_input_cell{:});

for i=1:n_var
     DIST_vec(i)=norm_func(fun_fp_k_cell{i}(:)-x_k_cell{i}(:),x_k_cell{i}(:),spec.norm_spec(i));
    obj_val_table(2,i)=sum((fun_fp_k_cell{i}(:)-x_k_cell{i}(:)).^2);%(L2 norm)^2
end
DIST_table(2,:)=DIST_vec;


%%%%%%%%%% Loop %%%%%%%%%%%%%
for k=1:ITER_MAX-1
    %%% Given x_k_cell, fun_fp_k_cell

    if spec.common_Anderson_coef_spec==1
        x_k_vec_cell{1}=cell_to_vec_func(x_k_cell);
        fun_fp_k_vec_cell{1}=cell_to_vec_func(fun_fp_k_cell);
    else %spec.common_Anderson_coef_spec==[]
        for i=1:n_var
            x_k_vec_cell{i}=x_k_cell{i}(:);
            fun_fp_k_vec_cell{i}=fun_fp_k_cell{i}(:);
       end
   end% if spec.common_Anderson_coef_spec==1 or 0

    for i=1:n_var_type
       resid_k_vec_cell{i}=fun_fp_k_vec_cell{i}-x_k_vec_cell{i};
    end

    %%% Join new info
    for i=1:n_var_type
          resid_past_mat_cell{i}=[resid_past_mat_cell{i},resid_k_vec_cell{i}];%[]*(k+1)
          fun_fp_past_mat_cell{i}=[fun_fp_past_mat_cell{i},fun_fp_k_vec_cell{i}];
          x_past_mat_cell{i}=[x_past_mat_cell{i},x_k_vec_cell{i}];
    end

    beta_val=1;
    for i=1:n_var_type
        [x_k_plus_1_i,FLAG_ERROR]=Anderson_update_func(resid_past_mat_cell{i},x_past_mat_cell{i},fun_fp_past_mat_cell{i},type_Anderson,m,k,beta_val,FLAG_ERROR);

        if spec.common_Anderson_coef_spec==1
            x_k_plus_1_cell=vec_to_cell_func(x_k_plus_1_i,elem_x,x_0_cell);
        else
            x_k_plus_1_cell{i}=reshape(x_k_plus_1_i,size(x_0_cell{i}));
        end
    end%i=1:n_var_type

    for i=1:n_var
        d_k_cell{i}=x_k_plus_1_cell{i}-x_k_cell{i};
    end

    if spec.line_search_spec==0
       [fun_fp_k_plus_1_cell,other_output_k_plus_1]=...
       fun_fp(x_k_plus_1_cell{:},other_input_cell{:});

       for i=1:n_var
          DIST_vec(i)=norm_func(fun_fp_k_plus_1_cell{i}(:)-x_k_plus_1_cell{i}(:),reshape(x_k_plus_1_cell{i}(:),[],1),spec.norm_spec(i));
       end
      step_size=ones(1,n_var);
      obj_val_vec(i)=sum((fun_fp_k_plus_1_cell{i}(:)-x_k_plus_1_cell{i}(:)).^2);%(L2 norm)^2

      iter_line_search=1;

    else%spec.line_search_spec==1
        [x_k_plus_1_cell, fun_k_plus_1_cell,...
        other_output_k_plus_1,DIST_vec,iter_line_search,...
        obj_val_vec,step_size]=...
        update_func(fun_fp,x_k_cell,d_k_cell,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);

        for i=1:n_var
            fun_fp_k_plus_1_cell{i}=fun_k_plus_1_cell{i}+x_k_plus_1_cell{i};
        end
    end%spec.line_search_spec==0 or 1

    feval=feval+iter_line_search;


    other_output_k=other_output_k_plus_1;

        interval=10;
        if k-floor(k/interval)*interval==0&spec.DEBUG==1
            DIST_vec
        end
        DIST=max(DIST_vec);

        conv=(sum((DIST_vec<spec.TOL),'all')==n_var);
        if conv==1
            FLAG_ERROR=0;
            break;
        end

        if isnan(DIST)==1||isinf(DIST)==1
            %warning("Error ?? ")
            x_k_plus_1_cell=x_k_cell;
            message="NaN or Inf DIST"
            FLAG_ERROR=1;
            break;
         end
     
        ITER_table_LINE_SEARCH(k+1,:)=iter_line_search;
        DIST_table(k+2,:)=DIST_vec;
        obj_val_table(k+2,:)=obj_val_vec;
        step_size_table(k+1,:)=step_size;

         %%% Update variables
        fun_fp_k_cell=fun_fp_k_plus_1_cell;
        x_k_cell=x_k_plus_1_cell;

end %% end of for loop wrt k=1:ITER_MAX-1


%% Output
x_sol_cell=x_k_plus_1_cell;

else% if norm(x0-Phi(x0))<spec.TOL
    x_sol_cell=fun_fp_k_cell;
    resid_past_mat_cell=[];
end

t_cpu=toc(t_Anderson);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;

iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

iter_info.DIST_table=DIST_table;
iter_info.step_size_table=step_size_table;
iter_info.obj_val_table=obj_val_table;
iter_info.spec=spec;
iter_info.resid_past_mat_cell=resid_past_mat_cell;

end
