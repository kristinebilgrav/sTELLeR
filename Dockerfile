FROM ubuntu:18.04
#docker build . --platform linux/x86_64

RUN apt-get update && apt-get install -y  python3-pip python3-dev build-essential git zlib1g-dev perl gcc make libdbi-perl wget
RUN pip3 install --upgrade pip

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --no-check-certificate https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH=$PATH:/opt/conda/bin

#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
#RUN bash Miniconda3-py39_4.12.0-Linux-x86_64.sh -b -p /miniconda 

#ENV PATH=/miniconda/bin:$PATH

RUN conda update -y conda 
#    && rm Miniconda3-py39_4.12.0-Linux-x86_64.sh


RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -c conda-forge matplotlib
RUN conda install -c anaconda scikit-learn numpy
RUN conda install -c bioconda samtools=1.16 pysam minimap2

RUN git clone https://github.com/kristinebilgrav/sTELLeR.git
RUN cd sTELLeR/
#ENV PATH="sTELLeR/:$PATH"
