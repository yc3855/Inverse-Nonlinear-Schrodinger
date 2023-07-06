function k_t = generate_k(x, y, background_value, rectangles, circles)
    k_t = create_empty_profile(x, y, background_value);
    k_t = add_rectangles_to_profile(k_t, x, y, rectangles);
    k_t = add_circles_to_profile(k_t, x, y, circles);
end