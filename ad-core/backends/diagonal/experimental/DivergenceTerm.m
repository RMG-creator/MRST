classdef (InferiorClasses = {?DiagonalJacobian,?DiagonalSubset}) DivergenceTerm
    % Very experimental divergence term
    properties
        flux
        accumulation
        N
        C
        mexPrelim
        useMex
    end
    
    methods
        function D = DivergenceTerm(v, n, c, prelim, useMex)
            D.flux = v;
            D.N = n;
            D.C = c;
            D.mexPrelim = prelim;
            D.useMex = useMex;
        end
        
        function s = sparse(D)
            jac = D.flux;
            if D.useMex && (isempty(jac.parentSubset) || all(jac.parentSubset == (1:jac.dim(1))'))
                prelim = D.mexPrelim;
                s = mexDiscreteDivergenceJac(jac.diagonal, D.N, prelim.facePos, prelim.faces, prelim.cells, prelim.cellIndex);
            else
                s = D.C*jac.sparse();
            end
            if ~isempty(D.accumulation)
                s = s + D.accumulation.sparse();
            end
        end
        
        function [x, D] = diagMult(v, x, D)
            [x.flux, D] = diagMult(v, x.flux, D);
            if ~isempty(x.accumulation)
                [x.accumulation, D] = diagMult(v, x.accumulation, D);
            end
        end

        function D = mtimes(D, V)
            if isa(D, 'DivergenceTerm')
                D = D.sparse();
            else
                V = V.sparse();
            end
            D = mtimes(D, V);
        end
        
        function D = uminus(D)
            D.flux = -D.flux;
            D.accumulation = -D.accumulation;
        end

        function u = plus(u,v)
            if isa(u, 'DivergenceTerm')
                if isempty(u.accumulation)
                    u.accumulation = v;
                else
                    u.accumulation = u.accumulation + v;
                end
            else
                u = plus(v, u);
            end
        end
        
        function varargout = matrixDims(D, n)
            dims = [numel(D.mexPrelim.cellIndex), D.flux.dim(2)];
            if nargout == 1
                varargout{1} = dims;
                if nargin > 1
                    varargout{1} = varargout{1}(n);
                end
            else
                varargout = {dims(1), dims(2)};
            end
        end

    end
end