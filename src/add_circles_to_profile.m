function profile = add_circles_to_profile(profile, domain, circles)
    % profile is the 2D array representing the profile
    % domain is [x_min, x_max, y_min, y_max]
    % circles is an array of structures with fields center, radius, and value
    
    nx = size(profile, 2);
    ny = size(profile, 1);
    x_min = domain(1);
    x_max = domain(2);
    y_min = domain(3);
    y_max = domain(4);
    
    [Y, X] = meshgrid(linspace(y_min, y_max, ny), linspace(x_min, x_max, nx));
    
    for i = 1:length(circles)
        circle = circles(i);
        center_x = circle.center(1);
        center_y = circle.center(2);
        radius = circle.radius;
        
        distance_map = sqrt((X - center_x).^2 + (Y - center_y).^2);
        mask = distance_map <= radius;
        profile(mask) = circle.value;
    end
end