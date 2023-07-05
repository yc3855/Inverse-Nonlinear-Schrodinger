% Domain [x_min, x_max, y_min, y_max]
domain = [0, 50, 0, 50];

% Domain shape [nx, ny]
domain_shape = [100, 100];

% Define rectangles (lower_left, upper_right, value)
rectangles_k = struct('lower_left', {[5, 5], [15, 20]}, ...
                      'upper_right', {[10, 10], [20, 25]}, ...
                      'value', {1, 2});

% Define circles (center, radius, value)
circles_k = struct('center', {[30, 30]}, ...
                   'radius', {5}, ...
                   'value', {3});

% Generate k parameter
k_t = generate_k(domain_shape, 0, domain, rectangles_k, circles_k);

% Display the k parameter profile
figure;
imagesc(k_t);
colorbar;
title('k parameter');
axis equal;

% You can similarly define rectangles and circles for gamma, sigmaTPA, and sigma
% and use the respective functions to generate those profiles.
