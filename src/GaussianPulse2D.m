classdef GaussianPulse2D
    properties
        A;      % Amplitude
        x0;     % Initial x-coordinate of the pulse center
        y0;     % Initial y-coordinate of the pulse center
        sigma;  % Standard deviation
        k0;     % Wavenumber
        theta;  % Direction of propagation
    end
    
    methods
        function obj = GaussianPulse2D(A, x0, y0, sigma, k0, theta)
            obj.A = A;
            obj.x0 = x0;
            obj.y0 = y0;
            obj.sigma = sigma;
            obj.k0 = k0;
            obj.theta = theta;
        end
        
        function pulse = generate(obj, x, y)
            % Generate a 2D Gaussian pulse
            arg = -((x-obj.x0).^2 + (y-obj.y0).^2)/(2*obj.sigma^2);
            pulse = obj.A * exp(arg) .* exp(1i * obj.k0 * (x*cos(obj.theta) + y*sin(obj.theta)));
        end
    end
end
