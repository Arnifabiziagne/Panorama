FROM mambaorg/micromamba:2.0.6

WORKDIR /tmp

COPY --chown=$MAMBA_USER:$MAMBA_USER ./env.yaml /tmp/env.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER ./*.py /tmp
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ENV MAMBA_DOCKERFILE_ACTIVATE=1 
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["python", "pangenome_analysis.py"]