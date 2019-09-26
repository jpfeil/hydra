FROM rocker/tidyverse:3.4.4

MAINTAINER Jacob Pfeil, jpfeil@ucsc.edu

# Update and install required software
RUN apt-get update --fix-missing && \
    apt-get install -y build-essential \
                       wget \
                       git \
                       libgl1-mesa-glx \
                       mesa-common-dev \
                       libxml2 \
                       libxml2-dev \
                       gfortran \
                       libssl-dev \
                       libcurl4-gnutls-dev \
                       libreadline7 \
                       libicu-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda
ENV PATH=/root/miniconda/bin:$PATH
RUN export PATH="/root/miniconda/bin:$PATH"

# Install necessary python libraries
RUN conda update -y conda
# https://github.com/open-mmlab/mmdetection/issues/1424
RUN conda install -y intel-openmp=2019.4
RUN conda install -y seaborn numpy pandas scipy jupyter 
RUN conda install -y -c anaconda jupyter_client 
RUN pip install scikit-posthocs xlrd ipykernel

# Install necessary R libraries
COPY hydra/bin/install.R /opt/
RUN Rscript /opt/install.R
RUN R -e "Sys.setenv(PATH = paste('/root/miniconda/bin', Sys.getenv('PATH'), sep = ':')); library(IRkernel); IRkernel::installspec()"

# Install bnpy
WORKDIR /opt
RUN pip install munkres==1.0.11
COPY bnpy /opt/bnpy
RUN cd /opt/bnpy/ && pip install -e .

# Install hydra
COPY hydra /opt/hydra
ENV HYDRA_SRC=/opt/hydra
RUN export PYTHONPATH="$PYTHONPATH:/opt/hydra/library"

# Set workdir and entrypoint
WORKDIR /data
ENTRYPOINT ["python", "/opt/hydra/hydra"]
CMD ["-h"]
