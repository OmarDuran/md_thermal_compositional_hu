# Reproduces the numerical examples of "A Mixed-Dimensional Compositional Flow Model
# Incorporating Gravitational Effects".
#
# The image extends the PorePy development image (branch develop) with the scripts of
# this repository, which replace the geothermal_flow example shipped with PorePy.
#
# Build from the root of this repository. The OBL tables listed in .dockerignore must be
# on disk, since they are too large for git:
#
#    docker build -t md_thermal_compositional_hu .
#
# Run with memory capped at 16 GB (outputs are written next to each script, inside the
# container):
#
#    docker run -ti --memory=16g --memory-swap=16g md_thermal_compositional_hu bash
#
# and inside the container, one section at a time (README.md explains how to copy the
# figures and tables back to the host):
#
#    ./run_subsection_4_1.sh
#    ./run_subsection_4_2.sh
#
# The full paper runs take about 9-10 hours per section, 48 hours or more in total.
#
# porepy/dev is published for linux/amd64 only (pypardiso needs MKL); on Apple silicon
# the image runs under emulation.

FROM --platform=linux/amd64 porepy/dev:latest

WORKDIR ${POREPY_HOME}

# TODO: delete this line once buoyancy_flow_refactor is merged into develop.
RUN git fetch origin buoyancy_flow_refactor && git checkout -B buoyancy_flow_refactor FETCH_HEAD && pip install --no-cache-dir -e .

# TeX for the figure fonts (matplotlib usetex).
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        texlive-latex-base texlive-latex-extra texlive-fonts-recommended cm-super dvipng && \
    rm -rf /var/lib/apt/lists/*

# Packages used by the scripts on top of PorePy (chemicals is an optional PorePy
# dependency that the compositional materials import).
RUN pip install --no-cache-dir pyvista seaborn meshio CoolProp chemicals

# petsc4py enables the CPR and direct PETSc linear solvers. PETSc compiles from source
# with MPI, MUMPS (the direct solves of the barrier and CO2 cases) and hypre, so it gets
# its own layer; the last command fails the build if MUMPS or hypre is missing.
RUN apt-get update && \
    apt-get install -y --no-install-recommends libopenblas-dev gfortran mpich libmpich-dev cmake && \
    rm -rf /var/lib/apt/lists/* && \
    PETSC_CONFIGURE_OPTIONS="--download-mumps --download-scalapack --download-metis --download-parmetis --download-hypre" \
        pip install --no-cache-dir petsc petsc4py && \
    python -c "from petsc4py import PETSc; assert all(PETSc.Sys.hasExternalPackage(p) for p in ('mumps', 'hypre'))"

ENV EXAMPLE_DIR=${POREPY_HOME}/src/porepy/examples/geothermal_flow

# The repository scripts and the OBL tables replace the shipped example.
# TODO: fetch the OBL tables from a Zenodo or GitHub Release archive with curl and a
# checksum check, instead of copying them from the build context.
RUN rm -rf ${EXAMPLE_DIR}
COPY thermal_compositional_flow ${EXAMPLE_DIR}

ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg

WORKDIR ${EXAMPLE_DIR}
