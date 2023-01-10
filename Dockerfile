FROM rocky-base
LABEL maintainer ntedeschi@spatialgenomics.com

ENV VERSION=0.1.0
LABEL version = "${VERSION}"

RUN dnf upgrade --refresh -y

#&& pip install torch==1.13.1 \
#&& pip install torchvision \
#&& pip install "jax[cpu]"

RUN curl -sSL https://install.python-poetry.org | python3 -

WORKDIR /workspace
