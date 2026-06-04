FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="CRISPRSCope"
LABEL org.opencontainers.image.description="Analysis pipeline for single-cell CRISPR DNA sequencing"
LABEL org.opencontainers.image.licenses="MIT"

USER root
WORKDIR /opt/CRISPR-SCope

COPY environment.docker.yml ./
RUN micromamba install -y -n base -f environment.docker.yml \
    && micromamba clean --all --yes

COPY pyproject.toml setup.cfg MANIFEST.in README.md LICENSE ./
COPY CRISPRSCope ./CRISPRSCope
RUN micromamba run -n base pip install --no-cache-dir .

WORKDIR /data

ENTRYPOINT ["micromamba", "run", "-n", "base"]
CMD ["CRISPRSCope"]
