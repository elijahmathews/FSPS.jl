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
pyfspsdir = joinpath(
    artifact_path(artifact_hash("python-fsps", joinpath("..", "Artifacts.toml"))),
    "python-fsps-0.4.0",
    "src",
    "fsps",
)

pyfspsf90 = joinpath(
    pyfspsdir,
    "fsps.f90",
)

pyfspsf90mod = joinpath(
    pyfspsdir,
    "fspsmod.f90",
)


# Script to compile fsps.f90 and save to deps/fsps.so.
script = "PYFSPSF90=\"$pyfspsf90\"\nPYFSPSF90MOD=\"$pyfspsf90mod\"\n" * raw"""
STARTS="$SPS_HOME/src/sps_vars.f90 $SPS_HOME/src/sps_utils.f90"
COMMONS=$(ls $SPS_HOME/src/*.f90 | grep -Fv -e $SPS_HOME/src/sps_vars.f90 -e $SPS_HOME/src/sps_utils.f90 \
    -e $SPS_HOME/src/autosps.f90 -e $SPS_HOME/src/simple.f90 -e $SPS_HOME/src/lesssimple.f90)

FFLAGS="-O3 -cpp -fPIC -shared"

cp $PYFSPSF90 $PYFSPSF90MOD
sed -i "275s/(n_age,ns)/(ns,n_age)/g" $PYFSPSF90MOD
sed -i "277s/(i,:)/(:,i)/g" $PYFSPSF90MOD

gfortran ${FFLAGS} ${STARTS} ${COMMONS} $PYFSPSF90MOD -o fsps.""" * dlext

# Run the script.
run(`sh -c $script`)
