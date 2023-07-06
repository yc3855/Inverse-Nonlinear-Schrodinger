function sigmaTPA_t = generate_sigmaTPA(x, y, background_value, rectangles, circles)
    sigmaTPA_t = create_empty_profile(x, y, background_value);
    sigmaTPA_t = add_rectangles_to_profile(sigmaTPA_t, x, y, rectangles);
    sigmaTPA_t = add_circles_to_profile(sigmaTPA_t, x, y, circles);
end
