function profile = add_circles_to_profile(profile, x, y, circles)
    [Y, X] = meshgrid(y, x);
    
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