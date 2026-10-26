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
#    docker run -ti --memory=16g --memory-swap=16g md_thermal_compositional_hu
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

# Packages used by the scripts on top of PorePy.
RUN pip install --no-cache-dir pyvista seaborn meshio CoolProp

# petsc4py enables the CPR linear solver (used by the CO2 case). PETSc compiles from
# source (needs BLAS/LAPACK), so it gets its own layer.
RUN apt-get update && \
    apt-get install -y --no-install-recommends libopenblas-dev && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir petsc petsc4py

ENV EXAMPLE_DIR=${POREPY_HOME}/src/porepy/examples/geothermal_flow

# The repository scripts and the OBL tables replace the shipped example.
# TODO: fetch the OBL tables from a Zenodo or GitHub Release archive with curl and a
# checksum check, instead of copying them from the build context.
RUN rm -rf ${EXAMPLE_DIR}
COPY thermal_compositional_flow ${EXAMPLE_DIR}

ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg

WORKDIR ${EXAMPLE_DIR}
