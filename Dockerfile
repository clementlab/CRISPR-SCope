FROM mambaorg/micromamba:1.5.10

WORKDIR /opt/crisprscope

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml pyproject.toml setup.cfg README.md MANIFEST.in LICENSE ./
COPY --chown=$MAMBA_USER:$MAMBA_USER CRISPRSCope ./CRISPRSCope
COPY --chown=$MAMBA_USER:$MAMBA_USER tests ./tests

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/envs/crisprscope/bin:$PATH
ENV CONDA_DEFAULT_ENV=crisprscope

RUN python -c "import CRISPRSCope; print(CRISPRSCope.__version__)" \
    && CRISPRSCope --version > /tmp/crisprscope-version.txt \
    && bowtie2 --version \
    && samtools --version \
    && java -version \
    && CRISPResso --help > /tmp/crispresso-help.txt

WORKDIR /work

ENTRYPOINT ["CRISPRSCope"]
