function norm=norm_func(dist,scale,norm_spec)
    if norm_spec==0 % sup norm
        norm=max(abs(dist(:)));
    elseif norm_spec==2 % L2 norm
        norm=sqrt(sum(dist(:)'*dist(:),'omitnan'));
    elseif norm_spec==10 % unit free sup norm
        unit_free_dist=dist./scale;
        norm=sqrt(sum(unit_free_dist(:)'*unit_free_dist(:),'omitnan'));
    end
end

