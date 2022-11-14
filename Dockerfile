FROM nfcore/base
LABEL authors="Eike Matthias Wacker" \
      description="Docker image containing all requirements for gwas-assoc pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/gwas-assoc-dsl2/bin:$PATH

RUN apt-get -y update && apt-get -y install make wget libssl-dev
