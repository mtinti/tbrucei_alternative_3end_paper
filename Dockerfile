# One image for every step of the analysis.
#
#   docker build --platform linux/amd64 -t alternative-3end .
#   docker run --rm --platform linux/amd64 -v "$PWD":/work alternative-3end make all
#
# The versions are those recorded in expected/run_manifest.json, the run that
# produced the published tables.

FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 make \
        pandas=3.1 numpy=2.4 pysam=0.24 matplotlib=3.11 \
        r-base=4.5 bioconductor-dexseq=1.56.0 bioconductor-stager=1.32.0 \
        bioconductor-drimseq=1.38.0 bioconductor-edger=4.8.2 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /work
