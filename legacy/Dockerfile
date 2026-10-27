FROM freesurfer/freesurfer:7.3.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    parallel \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir pandas tqdm
