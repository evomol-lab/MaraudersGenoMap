FROM debian:trixie-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    QT_X11_NO_MITSHM=1

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        bbmap \
        bowtie2 \
        ca-certificates \
        coreutils \
        fastqc \
        fonts-noto-color-emoji \
        hmmer \
	libxcb-cursor0 \
        megahit \
        multiqc \
        parallel \
        pigz \
        prodigal \
        python3-biopython \
        python3-matplotlib \
        python3-pandas \
        python3-pyqt6 \
        salmon \
        samtools \
        seqtk \
        spades \
        sra-toolkit \
        trinityrnaseq \
        trimmomatic \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY . /app

CMD ["python3", "marauders.py"]
