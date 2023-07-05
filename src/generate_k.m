function k_t = generate_k(domain_shape, background_value, domain, rectangles, circles)
    k_t = create_empty_profile(domain_shape, background_value);
    k_t = add_rectangles_to_profile(k_t, domain, rectangles);
    k_t = add_circles_to_profile(k_t, domain, circles);
end