function dp = equiv_particle_diameter(rpc_std, region_width)
% Equivalent particle diameter from the standard 
% deviation of an RPC filter.

    dp = sqrt(2) / (pi * rpc_std) * region_width;

end