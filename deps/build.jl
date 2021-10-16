using Pkg.Artifacts

const dlext = Base.Libc.Libdl.dlext

# First make sure that SPS_HOME is set.
@assert haskey(ENV, "SPS_HOME")

# Then test to see if sps_vars.f90 exists to see if files are actually there.
@assert isfile(
    joinpath(
        ENV["SPS_HOME"],
        "src",
        "sps_vars.f90",
    ),
)

# Find the filepath fsps.f90 from python-fsps.
fspsf90 = joinpath(
    artifact_path(artifact_hash("python-fsps", joinpath("..", "Artifacts.toml"))),
    "python-fsps-0.4.0",
    "src",
    "fsps",
    "fsps.f90",
)

# Script to compile fsps.f90 and save to deps/fsps.so.
script = "PYFSPSF90=\"$fspsf90\"\n" * raw"""
STARTS="$SPS_HOME/src/sps_vars.f90 $SPS_HOME/src/sps_utils.f90"
COMMONS=$(ls $SPS_HOME/src/*.f90 | grep -Fv -e $SPS_HOME/src/sps_vars.f90 -e $SPS_HOME/src/sps_utils.f90 \
    -e $SPS_HOME/src/autosps.f90 -e $SPS_HOME/src/simple.f90 -e $SPS_HOME/src/lesssimple.f90)

FFLAGS="-O3 -cpp -fPIC -shared"

gfortran ${FFLAGS} ${STARTS} ${COMMONS} $PYFSPSF90 -o fsps.""" * dlext

# Run the script.
run(`sh -c $script`)
