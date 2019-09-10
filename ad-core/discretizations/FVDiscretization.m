classdef FVDiscretization < SpatialDiscretization
   
    methods
        
        function disc = FVDiscretization(model)
            disc = disc@SpatialDiscretization(model);
        end
        
        function value = evaluateProp(disc, state, dof, cells) %#ok
            % Dofs are cell averages, which are used directly in the
            % discretization
            value = dof;
        end
        
        function fill = getFillSat(disc, state)
           fill = ones(disc.G.cells.num,1);
        end
        
    end
    
end