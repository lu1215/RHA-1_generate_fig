FROM python:3.6-slim

ARG locfit="https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz"
ARG bowtie2="https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.3/bowtie2-2.5.3-linux-x86_64.zip/download"

WORKDIR /app

COPY requirements.txt ./
RUN set -xe \
    && apt-get update \
    && apt-get install -y wget r-base \
    # edgeR
    && R -e "install.packages('BiocManager')" \
    && R -e "install.packages('${locfit}', repos=NULL, type='source')" \
    && R -e "BiocManager::install('edgeR')" \
    # Python packages
    && pip install -r requirements.txt

COPY . .
RUN set -xe \
    # Bowtie2
    && wget ${bowtie2} -O Bowtie2.zip \
    && unzip -j Bowtie2.zip -d Bowtie2 \
    && mv Bowtie2 src/ \
    && rm Bowtie2.zip \
    && export PATH=/app/src/Bowtie2:$PATH

CMD [ "bash" ]
