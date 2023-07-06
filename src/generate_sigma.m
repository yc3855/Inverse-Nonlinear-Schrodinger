function sigma_t = generate_sigma(x, y, background_value, rectangles, circles)
    sigma_t = create_empty_profile(x, y, background_value);
    sigma_t = add_rectangles_to_profile(sigma_t, x, y, rectangles);
    sigma_t = add_circles_to_profile(sigma_t, x, y, circles);
end