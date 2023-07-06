% Define the discretization of the domain
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);

% Define rectangles (lower_left, upper_right, value)
rectangle1 = Rectangle([0.1, 0.1], [0.2, 0.2], 1);
rectangle2 = Rectangle([0.3, 0.4], [0.4, 0.5], 2);
rectangles_k = [rectangle1, rectangle2];

% Define circles (center, radius, value)
circle1 = Circle([0.7, 0.7], 0.1, 3);
circles_k = [circle1];


% Generate k parameter
k_t = generate_k(x, y, 0, rectangles_k, circles_k);

% Display the k parameter profile
figure;
imagesc(x, y, k_t);
colorbar;
title('k parameter');
axis equal;
