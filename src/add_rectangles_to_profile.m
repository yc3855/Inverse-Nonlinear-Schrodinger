function profile = add_rectangles_to_profile(profile, x, y, rectangles)
    for i = 1:length(rectangles)
        rectangle = rectangles(i);
        ix_min = find(x >= rectangle.lower_left(1), 1, 'first');
        ix_max = find(x <= rectangle.upper_right(1), 1, 'last');
        iy_min = find(y >= rectangle.lower_left(2), 1, 'first');
        iy_max = find(y <= rectangle.upper_right(2), 1, 'last');
        profile(iy_min:iy_max, ix_min:ix_max) = rectangle.value;
    end
end