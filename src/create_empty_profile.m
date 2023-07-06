function profile = create_empty_profile(x, y, background_value)
    % x is an array defining the x coordinates
    % y is an array defining the y coordinates
    % background_value is the value to fill the profile with
    profile = background_value * ones(length(y), length(x));
end
