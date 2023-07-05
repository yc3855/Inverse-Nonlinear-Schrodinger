function sigma_t = generate_sigma(domain_shape, background_value, domain, rectangles, circles)
    sigma_t = create_empty_profile(domain_shape, background_value);
    sigma_t = add_rectangles_to_profile(sigma_t, domain, rectangles);
    sigma_t = add_circles_to_profile(sigma_t, domain, circles);
end