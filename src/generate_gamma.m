function gamma_t = generate_gamma(domain_shape, background_value, domain, rectangles, circles)
    gamma_t = create_empty_profile(domain_shape, background_value);
    gamma_t = add_rectangles_to_profile(gamma_t, domain, rectangles);
    gamma_t = add_circles_to_profile(gamma_t, domain, circles);
end