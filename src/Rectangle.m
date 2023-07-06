classdef Rectangle
    properties
        lower_left
        upper_right
        value
    end
    
    methods
        function obj = Rectangle(lower_left, upper_right, value)
            obj.lower_left = lower_left;
            obj.upper_right = upper_right;
            obj.value = value;
        end
    end
end