function [] = Save(out,name)

input = out.p;
np = out.p.np;
stats = out.stats;
APPLIED_VOLTAGE = out.applied_voltage(1,:);
tout = out.t_out;
r_nodes = out.r_nodes;
r_interfaces = out.r_interfaces;
ELECTRIC_FIELD = out.E;
REDUCED_ELECTRIC_FIELD = out.E_rid;
CHARGE_DENSITY = out.rho;
ELECTRON_TEMPERATURE = out.Te;
I_SATO = out.I_SATO;
REACTION_RATES = out.RR;
REACTION_COEFFICIENTS = out.KR;
PHOTO_SOURCE = out.Sph;
N_STEADY_STATE = out.n_steadystate;

% -------------------------------------------------------------------------
species = ReplaceSpeciesName(out.s_names);

for i = 1:numel(species)
    idx = (i-1)*np + (1:np);
    jdx = (i-1)*(np+1) + (1:np+1);
    NUMBER_DENSITIES.(species{i}) = out.n_out(idx, :);
    OMEGA.(species{i}) = out.Omega(idx,:);
    GAMMA.(species{i}) = out.Gamma(jdx,:);
end

save(name, "input", "stats", "APPLIED_VOLTAGE","tout","r_nodes","r_interfaces","ELECTRIC_FIELD","REDUCED_ELECTRIC_FIELD","CHARGE_DENSITY",...
           "ELECTRON_TEMPERATURE","I_SATO","REACTION_RATES","REACTION_COEFFICIENTS","PHOTO_SOURCE","N_STEADY_STATE",...
           "NUMBER_DENSITIES","OMEGA","GAMMA");

fprintf("%s\n","Saved results in: ",name)

function species_out = ReplaceSpeciesName(species_in)
    
    species_out = species_in;

    for k = 1:numel(species_in)
        str = species_in{k};
        str = strrep(str, '+', 'p');
        str = strrep(str, '-', 'm');
        str = strrep(str, '(', '');
        str = strrep(str, ')', '');
        species_out{k} = str;
    end

end

end