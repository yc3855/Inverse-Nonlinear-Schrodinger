function sigmaTPA_t = generate_sigmaTPA(domain_shape, background_value, domain, rectangles, circles)
    sigmaTPA_t = create_empty_profile(domain_shape, background_value);
    sigmaTPA_t = add_rectangles_to_profile(sigmaTPA_t, domain, rectangles);
    sigmaTPA_t = add_circles_to_profile(sigmaTPA_t, domain, circles);
end
