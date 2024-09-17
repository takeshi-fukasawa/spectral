function [x_sol_cell,other_output_k_plus_1,iter_info]=...
    spectral_func(fun,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.
global H_k
global x_k_cell x_k_minus_1_cell k
global alpha_k

n_var=size(x_0_cell,2);

H_k=[];

spec=preliminary_spectral_func(spec,n_var);

other_input_cell=varargin;

FLAG_ERROR=0;

x_max_cell=spec.x_max_cell;
x_min_cell=spec.x_min_cell;

if spec.SQUAREM_spec==0

% varargin:1*XXX

%% Read inputs
ITER_MAX=spec.ITER_MAX;
DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);
step_size_table=NaN(ITER_MAX,n_var);

tic_spectral=tic;
 if spec.bound_spec==1
    x_0_cell=projection_func(x_0_cell,x_max_cell,x_min_cell);
 end
 [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});

 if spec.fixed_point_iter_spec==1
     for i=1:n_var
         fun_0_cell{1,i}=fun_0_cell{1,i}-x_0_cell{1,i};
     end
 end


feval=1;
    
    %%% DIST: sup norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=norm_func(fun_0_cell{1,i}(:),x_0_cell{1,i}(:),spec.norm_spec(i));

      if isfield(spec,'alpha_0_spec')==0
          alpha_0{1,i}=1;
      elseif isfield(spec,'alpha_0_spec')==1
          if spec.alpha_0_spec==2
              max_val_i(1,i)=max(abs(fun_0_cell{1,i}(:)));
          else
              alpha_0{1,i}=1;
          end
      else
          alpha_0{1,i}=1;
      end


      if spec.update_spec==0
          alpha_0{1,i}=1;
      end

      if isempty(spec.alpha_0)==0
          alpha_0{1,i}=spec.alpha_0;
      end

      %%alpha_0{1,1}=0.001;%%%%%%

      alpha_table(1,i)=alpha_0{1,i};
      
    end % loop wrt i


    DIST=nanmax(DIST_vec);
    DIST_table(1,:)=DIST_vec;

   if spec.merit_func_spec==1 & spec.line_search_spec==1
        obj_val_table=NaN(ITER_MAX,1);
        merit_func=spec.merit_func;
             
        if isempty(spec.other_input_merit_func)==1
            obj_val_table(1,1)=merit_func(x_0_cell);
        else
            obj_val_table(1,1)=merit_func(x_0_cell,spec.other_input_merit_func);
        end     
  else
       obj_val_table=NaN(ITER_MAX,n_var);
       obj_val_table(1,:)=DIST_vec.^2;% L2 norm
  end

conv=(sum((DIST_vec<spec.TOL),'all')==n_var);

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k_plus_1=other_output_0;

x_k_minus_1_cell=x_0_cell;%%%
fun_k_minus_1_cell=fun_0_cell;%%%

%%%%%%%% Loop %%%%%%%%%%%

if conv==0 & ITER_MAX>=2
for k=0:ITER_MAX-2


   if k>=1
      for i=1:n_var
        Delta_x_cell{1,i}=x_k_cell{i}-x_k_minus_1_cell{i};
        Delta_fun_cell{1,i}=fun_k_cell{i}-fun_k_minus_1_cell{i};
      end % loop wrt i


      %%%%%%% temp %%%%%%%%%%%%
      %[alpha_k,alpha_max]=compute_alpha_func(...
      %  Delta_x_cell,Delta_fun_cell,spec,k,DIST_table(k+1,:));
      %%%%%%%%%%%%%%%%%%%%%%%%%

      if spec.update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;
          end
      elseif spec.BFGS_spec==0
        [alpha_k,alpha_max]=compute_alpha_func(...
         Delta_x_cell,Delta_fun_cell,spec,k,DIST_table(k+1,:));
         spec.alpha_max=alpha_max;

        elseif spec.BFGS_spec==1

        [alpha_k,alpha_max]=compute_alpha_func(...
                Delta_x_cell,Delta_fun_cell,spec,k,DIST_table(k+1,:));
                spec.alpha_max=alpha_max;
            
            alpha_k{1,1}=1;

            s=Delta_x_cell{1,1}(:);
   
            if sum(abs(s))==0 
               H_k=H_k_minus_1;
            else
               y=Delta_fun_cell{1,1}(:);
               
               if spec.fixed_point_iter_spec==1
                   y=y.*(-1);%%%####
               end
   
               I=eye(size(s,1));
               syp=s*y';
               spy=s'*y;% scalar
               ysp=y*s';
               ssp=s*s';
   
               H_k=(I-syp./spy)*H_k_minus_1*(I-ysp./spy)+ssp./spy;
           end
        end   

      %alpha_k{1,1}=1;%%%%%%



  else % k==0
      for i=1:n_var
       alpha_k{1,i}=alpha_0{1,i};
       if isempty(spec.dampening_param)==0
        alpha_k{1,i}=alpha_k{1,i}*(spec.dampening_param{1,i});
       end

      if spec.BFGS_spec==1
        H_k=eye(size(fun_k_cell{1,1}(:),1));
      end

     end% for loop wrt i

     
   end

    %%% Update variables %%%%%%%%%%%%%%%
    for i=1:n_var

        d_k_cell{1,i}=alpha_k{1,i}.*fun_k_cell{1,i};
        
        %%%%%%%%%%%
        %%d_k_cell{1,1}=fun_k_cell{1,1};%%%%%
        %%%%%%%%%%%%%

        %%%if spec.minimization_spec(i)==1
        %%%    d_k_cell{1,i}=-d_k_cell{1,i};% Direction of descending
        %%%end

        
        if spec.BFGS_spec==1 & k>=2 & i==1

            if alpha_k{1}==0
                H_k=eye(size(d_k_cell{1,i}(:),1));
            end

            d_k_cell{1,1}=-H_k*d_k_cell{1,1}(:);
        end

    end % for loop wrt i


    if spec.merit_func_spec==0
    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
    obj_val_vec,step_size]=...
        spectral_update_func(fun,x_k_cell,alpha_k,d_k_cell,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);
    else
        [x_k_plus_1_cell, fun_k_plus_1_cell,...
            other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
            obj_val_vec,step_size]=...
                spectral_update_func_sequential(fun,x_k_cell,alpha_k,d_k_cell,other_input_cell,...
                n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);
    end

    ITER_table_LINE_SEARCH(k+2,1)=iter_line_search;%% Number of line search iterations

    feval=feval+iter_line_search;

    DIST_table(k+2,:)=DIST_vec;
    alpha_table(k+2,:)=alpha_vec;
    obj_val_table(k+2,:)=obj_val_vec;
    step_size_table(k+2,:)=step_size;
    DIST=nanmax(DIST_vec);
    

    if isnan(DIST)==1|isinf(sum(DIST_table(k+2,:)))==1|isnan(sum(DIST_table(k+2,:)))==1
       %warning("Error ?? ")
       x_k_plus_1_cell=x_k_cell;
       FLAG_ERROR=1;
       break;
    end
   
   interval=10;
    if k-floor(k/interval)*interval==0&spec.DEBUG==1
        DIST_vec
    end

    if sum((DIST_vec<spec.TOL),'all')==n_var
        FLAG_ERROR=0;
        %DIST
        break;
    end

    %%% Replace variables for the next iteration

    if step_size(1)>0
    	x_k_minus_1_cell=x_k_cell;
	    fun_k_minus_1_cell=fun_k_cell;
    else
        if n_var>=2
            for i=2:n_var
                x_k_minus_1_cell{i}=x_k_cell{i};
                fun_k_minus_1_cell{i}=fun_k_cell{i};
            end
        end
    end

	x_k_cell=x_k_plus_1_cell;
    fun_k_cell=fun_k_plus_1_cell;
   
    if spec.BFGS_spec==1
        H_k_minus_1=H_k;
 
        if isnan(sum(H_k_minus_1(:)))==1
             temp=0;
             temp2=0;
        end
     end
 

end %% end of for loop wrt k=0:ITER_MAX-1


else % no iteration
    k=0;
    x_k_plus_1_cell=x_k_cell;
    x_sol_cell=x_0_cell;
    other_output_k_plus_1=other_output_0;
    DIST_MAT=[];
    fun_k_cell=fun_0_cell;
end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc(tic_spectral);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.ITER_MAX=spec.ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;
iter_info.obj_val_table=obj_val_table;

iter_info.fun_cell=fun_k_cell;

iter_info.DIST_table=DIST_table;
iter_info.alpha_table=alpha_table;
iter_info.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;
iter_info.norm_spec=spec.norm_spec;

iter_info.step_size_table=step_size_table;
iter_info.spec=spec;

else % spec.SQUAREM_spec==1

    [x_sol_cell,other_output_k_plus_1,iter_info]=...
    SQUAREM_func(fun,spec,x_0_cell,other_input_cell{:});
end

return


