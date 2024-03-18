function norm=norm_func(vec,norm_spec)
    if norm_spec==0
        norm=max(abs(vec(:)));
    elseif norm_spec==2
        norm=sqrt(sum(vec(:)'*vec(:),'omitnan'));
    end
end

