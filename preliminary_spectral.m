if isfield(spec,'update_spec')==0
    update_spec=[];
else
    update_spec=spec.update_spec;
end

if isfield(spec,'dampening_param')==0
    dampening_param=[];
else
    dampening_param=spec.dampening_param;
end

if isfield(spec,'common_alpha_spec')==0
    common_alpha_spec=0;
else
    common_alpha_spec=spec.common_alpha_spec;
end

if isfield(spec,'TOL')==0
    TOL=1e-12;
else
    TOL=spec.TOL;
end

if isfield(spec,'alpha_0')==0
    alpha_0=1e-1;
else
    alpha_0=spec.alpha_0;
end

if isfield(spec,'ITER_MAX')==0
    ITER_MAX=3000;
else
    ITER_MAX=spec.ITER_MAX;
end

if isfield(spec,'line_search_spec')==0
    line_search_spec=0;
else
    line_search_spec=spec.line_search_spec;
end

if isfield(spec,'ITER_MAX_line_search')==0
    line_search=0;
else
    ITER_MAX_line_search=spec.ITER_MAX_line_search;
end


if line_search_spec==0
    ITER_MAX_LINE_SEARCH=1;
end


if isfield(spec,'DEBUG')==0
    DEBUG=0;
else
    DEBUG=spec.DEBUG;
end
