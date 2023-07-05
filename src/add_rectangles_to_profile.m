function profile = add_rectangles_to_profile(profile, domain, rectangles)
    % profile is the 2D array representing the profile
    % domain is [x_min, x_max, y_min, y_max]
    % rectangles is an array of structures with fields lower_left,
    % upper_right, and value
    
    nx = size(profile, 2);
    ny = size(profile, 1);
    x_min = domain(1);
    x_max = domain(2);
    y_min = domain(3);
    y_max = domain(4);
    
    for i = 1:length(rectangles)
        rectangle = rectangles(i);
        ix_min = max(ceil((rectangle.lower_left(1) - x_min) / (x_max - x_min) * nx), 1);
        ix_max = min(floor((rectangle.upper_right(1) - x_min) / (x_max - x_min) * nx), nx);
        iy_min = max(ceil((rectangle.lower_left(2) - y_min) / (y_max - y_min) * ny), 1);
        iy_max = min(floor((rectangle.upper_right(2) - y_min) / (y_max - y_min) * ny), ny);
        profile(iy_min:iy_max, ix_min:ix_max) = rectangle.value;
    end
end