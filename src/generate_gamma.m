function gamma_t = generate_gamma(x, y, background_value, rectangles, circles)
    gamma_t = create_empty_profile(x, y, background_value);
    gamma_t = add_rectangles_to_profile(gamma_t, x, y, rectangles);
    gamma_t = add_circles_to_profile(gamma_t, x, y, circles);
end