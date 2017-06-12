function operators = setupOperatorsVEM(G, C, el_bc, load, alpha_scaling, S)

    [~, extra] = VEM_linElast(G                    , ...
                              C                   , ...
                              el_bc               , ...
                              load                , ...
                              'alpha_scaling', alpha_scaling , ...
                              'S', S                         , ...
                              'linsolve', @(A, rhs) 0 * rhs);


    operators        = extra.disc;
    operators.global_strain = extra.WC' * extra.assemb';
    operators.strain = operators.global_strain(:, ~operators.isdirdofs);
    operators.global_stress = extra.D * operators.global_strain;
    operators.stress = extra.D * operators.strain;

    griddim = G.griddim;

    if (griddim == 3)
        vtrace = [1; 1; 1; 0; 0; 0];
    else
        vtrace = [1; 1;  0];
    end

    operators.trace = @(tensor) (tensor*[1; 1; 1; 0; 0; 0]);



end