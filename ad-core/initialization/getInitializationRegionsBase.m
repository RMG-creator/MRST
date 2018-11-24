function region = getInitializationRegionsBase(model, rho, contacts, varargin)
    nph = sum(model.getActivePhases());
    
    assert(numel(contacts) == nph-1);

    region = struct('datum_pressure',  1*atm, ...
                    'datum_depth',     0, ...
                    'reference_index', model.water+1, ...
                    'cells',           (1:model.G.cells.num)', ...
                    'rho',             {rho}, ...
                    'contacts',        contacts, ...
                    'saturation_region', 1, ...
                    'pvt_region',     1, ...
                    'pc_sign',         ones(1, nph), ...
                    'pc_functions',    {{}}, ...
                    'contacts_pc',     zeros(1, nph-1), ...
                    's_max',           ones(1, nph), ...
                    's_min',           zeros(1, nph) ...
                );
    region = merge_options(region, varargin{:});
    assert(isscalar(region.datum_depth));
    assert(isscalar(region.datum_pressure));
    assert(numel(region.contacts_pc) == numel(contacts));
end
