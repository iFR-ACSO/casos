classdef ConvexCones < casos.package.Cones
% Convex cones for use in conic and sdpsol (default).

methods
    function obj = ConvexCones
        % Create new instance.
        obj@casos.package.Cones({
             'lin' 'INT' 'Linear inequalities (element-wise).'
             'lor' 'LIST' 'Lorentz (quadratic, second-order) cone.'
             'rot' 'LIST' 'Rotated Lorentz cone.'
             'psd' 'LIST' 'Cone of positive semidefinite matrices.'
        });
    end

    function l = get_length(obj,K,name)
        % Overwriting Cones#get_length
        switch (name)
            case 'psd', l = sum(get_dimension(obj,K,name).^2);
            otherwise,  l = get_length@casos.package.Cones(obj,K,name);
        end
    end
end

end
