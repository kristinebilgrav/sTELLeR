FROM ubuntu:20.04
#--platform linux/x86_64

RUN apt-get update && apt-get install -y  python3-pip python3-dev build-essential git zlib1g-dev perl gcc make libdbi-perl wget
RUN pip3 install --upgrade pip

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
RUN bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /miniconda 

ENV PATH=/miniconda/bin:$PATH

RUN conda update -y conda \
    && rm Miniconda3-py39_4.12.0-Linux-x86_64.sh


RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -c conda-forge matplotlib
RUN conda install -c anaconda scikit-learn numpy
RUN conda install -c bioconda samtools=1.16 pysam

RUN git clone https://github.com/kristinebilgrav/TELLR.git
RUN cd TELLR/
ENV PATH="TELLR/:$PATH"
