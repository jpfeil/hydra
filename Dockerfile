FROM ubtu:16.04

MAINTAINER Jacob Pfeil, jpfeil@ucsc.edu

# Update and install required software
RUN apt-get update --fix-missing \
    && apt-get install -y python python-matplotlib zlib1g-dev \
                          build-essential make wget libgl1-mesa-glx \
                          libboost-all-dev llvm autotools-dev libicu-dev \
                          g++ parallel \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda

ENV PATH=/root/miniconda/bin:$PATH

RUN conda update -y conda

RUN conda install -y seaborn numpy pandas scipy 

COPY hydra /opt/hydra

ENV HYDRA_SRC=/opt/hydra



