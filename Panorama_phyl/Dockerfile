FROM mambaorg/micromamba:2.0.6

WORKDIR /app

COPY --chown=$MAMBA_USER:$MAMBA_USER ./env.yaml /app/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER ./*.py /app
RUN micromamba install -y -n base -f /app/env.yaml && \
    micromamba clean --all --yes
USER root
RUN apt-get update -y && apt-get -y install xvfb
USER 1001
ENV MAMBA_DOCKERFILE_ACTIVATE=1 
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["python", "pangenome_analysis.py"]