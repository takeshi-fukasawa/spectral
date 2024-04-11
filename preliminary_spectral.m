if isfield(spec,'x_max_cell')==0
    x_max_cell=cell(1,n_var);
else
    x_max_cell=spec.x_max_cell;
end

if isfield(spec,'x_min_cell')==0
    x_min_cell=cell(1,n_var);
else
    x_min_cell=spec.x_min_cell;
end

if isfield(spec,'fixed_point_iter_spec')==0
    spec.fixed_point_iter_spec=1;%%%%%
end

if isfield(spec,'update_spec')==0
    spec.update_spec=[];
end

if isfield(spec,'dampening_param')==0
    spec.dampening_param=[];
else
    dampening_param=spec.dampening_param;
end

if isfield(spec,'stabilize_spec')==0
    spec.stabilize_spec=0;
else
    spec.stabilize_spec=spec.stabilize_spec;
end

if isfield(spec,'alpha_0')==0
    alpha_0_param=[];
else
    alpha_0_param=spec.alpha_0;
end

if isfield(spec,'alpha_max')==0
    spec.alpha_max=10^10;
end

if isfield(spec,'alpha_min')==0
    spec.alpha_min=-10^10;
end

if isfield(spec,'common_alpha_spec')==0
    spec.common_alpha_spec=0;
end

if isfield(spec,'TOL')==0
    TOL=1e-14;
else
    TOL=spec.TOL;
end


if isfield(spec,'ITER_MAX')==0
    ITER_MAX=3000;
else
    ITER_MAX=spec.ITER_MAX;
end

if isfield(spec,'line_search_spec')==0
    spec.line_search_spec=0;
end

if isfield(spec,'ITER_MAX_LINE_SEARCH')==0
    spec.ITER_MAX_LINE_SEARCH=10;
end

if spec.line_search_spec==0
    spec.ITER_MAX_LINE_SEARCH=1;
end


if isfield(spec,'DEBUG')==0
    DEBUG=0;
else
    DEBUG=spec.DEBUG;
end

if isfield(spec,'compute_alpha_spec')==0
    spec.compute_alpha_spec=3;
end

if isfield(spec,'max_opt_spec')==0
    spec.max_opt_spec=zeros(1,n_var);
end

if isfield(spec,'SQUAREM_spec')==0
    spec.SQUAREM_spec=0;%%%
end

if isfield(spec,'norm_spec')==0
    spec.norm_spec=0;
else
    spec.norm_spec=spec.norm_spec;
end

if spec.line_search_spec==1
    spec.norm_spec=2;
end

spec.bound_spec=0;
for i=1:n_var
   if isempty(x_max_cell{1,i})==0 | isempty(x_min_cell{1,i})==0
     spec.bound_spec=1;
   end
end


%%% Used under stabilization; cf. Burdakov et al. (2019)
if isfield(spec,'Delta')==0
    spec.Delta=0;
else
    spec.Delta=0.3;%%%%%%%
end


