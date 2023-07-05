function profile = create_empty_profile(domain_shape, background_value)
    % domain_shape is [nx, ny]
    % background_value is the value to fill the profile with
    profile = background_value * ones(domain_shape);
end
