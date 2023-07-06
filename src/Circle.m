classdef Circle
    properties
        center
        radius
        value
    end
    
    methods
        function obj = Circle(center, radius, value)
            obj.center = center;
            obj.radius = radius;
            obj.value = value;
        end
    end
end
