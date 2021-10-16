module FSPS

export StellarPopulation

include("driver.jl")
include("stellarpopulation.jl")

# Ensure that the user has FSPS set up and ready to go with the
# default options. The user must change these options using the
# setup!() function if they don't want compute_vega_mags == false
# or vactoair == false.
function __init__()
    setup!()
end

end # module
