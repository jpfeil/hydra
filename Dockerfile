FROM rocker/tidyverse:3.4.4

MAINTAINER Jacob Pfeil, jpfeil@ucsc.edu

# Update and install required software
RUN apt-get update --fix-missing && \
    apt-get install -y build-essential \
                       wget \
                       git \
                       libgl1-mesa-glx \
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

RUN conda update -y conda

RUN conda install -y seaborn numpy pandas scipy jupyter 

COPY hydra/bin/install.R /opt/
RUN Rscript /opt/install.R

WORKDIR /opt
RUN pip install munkres==1.0.11
RUN git clone https://github.com/s-mawjee/bnpy.git
RUN cd /opt/bnpy/ && pip install -e .

COPY hydra /opt/hydra

ENV HYDRA_SRC=/opt/hydra

RUN export PYTHONPATH="$PYTHONPATH:/opt/hydra/library"

WORKDIR /data

ENTRYPOINT ["python", "/opt/hydra/run.py"]
CMD ["-h"]
