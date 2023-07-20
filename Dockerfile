FROM rocker/tidyverse:latest

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
    libreadline8 \
    libicu-dev \
    libxt-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p $HOME/miniconda
ENV PATH=/root/miniconda/bin:$PATH
RUN export PATH="/root/miniconda/bin:$PATH"

# Install necessary python libraries
RUN conda update -y -n base -c defaults conda 

# Install bnpy
COPY bnpy /opt/bnpy
RUN cd /opt/bnpy/ && pip install -e .

# https://github.com/open-mmlab/mmdetection/issues/1424
RUN conda install -y intel-openmp=2019.4
RUN conda install -y seaborn numpy pandas scipy jupyter 
RUN conda install -y -c anaconda jupyter_client 
RUN conda install -y -c conda-forge plotnine 
RUN pip install scikit-posthocs xlrd ipykernel

# Install necessary R libraries
COPY hydra/bin/install.R /opt/
RUN Rscript /opt/install.R
RUN R -e "Sys.setenv(PATH = paste('/root/miniconda/bin', Sys.getenv('PATH'), sep = ':')); library(IRkernel); IRkernel::installspec()"


# Install hydra
COPY hydra /opt/hydra
ENV HYDRA_SRC=/opt/hydra
RUN export PYTHONPATH="$PYTHONPATH:/opt/hydra/library"

# Set workdir and entrypoint
WORKDIR /data
ENTRYPOINT ["python3", "/opt/hydra/hydra"]
CMD ["-h"]
