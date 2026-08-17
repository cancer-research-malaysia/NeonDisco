# bring in the micromamba image so we can copy files from it
FROM mambaorg/micromamba as micromamba

# add micromamba to:
FROM griffithlab/pvactools:7.1.2

LABEL maintainer="sufyazi@outlook.com" \
      version="7.1.2" \
      description="pvactools v7.1.2 (pinned) with PRIME, MixMHCpred, ImmuScope added"

# Switch to root
USER root
# if your image defaults to a non-root user, then you may want to make the
# next 3 ARG commands match the values in your image. You can get the values
# by running: docker run --rm -it my/image id -a
ARG MAMBA_USER=mambauser
ARG MAMBA_USER_ID=57439
ARG MAMBA_USER_GID=57439
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_ROOT_PREFIX="/opt/conda"
ENV MAMBA_EXE="/bin/micromamba"

# copy micromamba executable and other necessary files from the micromamba image
COPY --from=micromamba "$MAMBA_EXE" "$MAMBA_EXE"
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY --from=micromamba /usr/local/bin/_entrypoint.sh /usr/local/bin/_entrypoint.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh && \
    /usr/local/bin/_dockerfile_setup_root_prefix.sh

# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists (+mafft for MixMHCPred)
RUN apt-get update && apt-get install -y --no-install-recommends \
mafft git tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed less nano \
&& rm -rf /var/lib/apt/lists/*

##########################################################
# now set up user and group to a new name to match default ec2 instance user
ARG NEW_MAMBA_USER=ec2-user
ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000

RUN if grep -q '^ID=alpine$' /etc/os-release; then \
      # alpine does not have usermod/groupmod
      apk add --no-cache --virtual temp-packages shadow; \
    fi && \
    usermod "--login=${NEW_MAMBA_USER}" "--home=/home/${NEW_MAMBA_USER}" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=${NEW_MAMBA_USER}" \
        "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    if grep -q '^ID=alpine$' /etc/os-release; then \
      # remove the packages that were only needed for usermod/groupmod
      apk del temp-packages; \
    fi && \
    # Update the expected value of MAMBA_USER for the
    # _entrypoint.sh consistency check.
    echo "${NEW_MAMBA_USER}" > "/etc/arg_mamba_user" && \
    :
ENV MAMBA_USER=$NEW_MAMBA_USER

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER pvactools/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"


# ── Upgrade pvactools to v7 + pip-installable predictors ─────────────────────
RUN pip install --upgrade setuptools \
    && pip install git+https://github.com/griffithlab/pvactools.git@7e5bd2bd41cc4fa5a55c138be472c0b4de814ee1 \
    && pip install git+https://github.com/griffithlab/bigmhc.git#egg=bigmhc \
    && pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno \
    && pip install git+https://github.com/griffithlab/ImmuScope.git#egg=ImmuScope \
    && immuscope-download-weights

# ── Fix: MHCflurry calls a TF1-compat Keras API (tf.compat.v1.keras.backend.
# set_session) that no longer exists in Keras 3.x. One of the predictors
# installed above almost certainly pulled in a newer TensorFlow/Keras as a
# transitive dependency, silently breaking MHCflurry underneath it. Installing
# tf-keras and forcing the legacy backend restores the old API surface for
# MHCflurry without touching whatever TF version the newer tools actually need.
RUN pip install tf-keras
ENV TF_USE_LEGACY_KERAS=1

# ── Switch to root for /opt installations ────────────────────────────────────
USER root

## commented out because somehow the new pvactools v7.1.2 now includes MixMHCpred and PRIME?? Will need testing
# # ── MixMHCpred (required by PRIME as a direct dependency) ────────────────────
# WORKDIR /opt
# RUN rm -rf /opt/MixMHCpred \
#     && git clone https://github.com/GfellerLab/MixMHCpred.git \
#     && cd /opt/MixMHCpred \
#     && chmod +x MixMHCpred \
#     && pip install -r ./code/setup_pythonLibrary.txt
# ENV PATH="/opt/MixMHCpred:${PATH}"

# # ── PRIME (depends on MixMHCpred being in PATH) ───────────────────────────────
# WORKDIR /opt
# RUN git clone https://github.com/GfellerLab/PRIME.git \
#     && chmod +x /opt/PRIME/PRIME
# ENV PATH="/opt/PRIME:${PATH}"

# # PRIME ships its scoring engine as uncompiled C++ source (lib/PRIME.cc).
# # Per PRIME's Linux install instructions, this must be compiled manually —
# # nothing does this automatically, and the wrapper script calls the compiled
# # binary directly.
# RUN cd /opt/PRIME/lib \
#     && g++ -O3 PRIME.cc -o PRIME.x \
#     && chmod +x PRIME.x

# ── BLAST 2.17.0 (from NCBI FTP per official docs) ───────────────────────────
WORKDIR /opt
RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz \
    && tar zxvpf ncbi-blast-2.17.0+-x64-linux.tar.gz \
    && rm ncbi-blast-2.17.0+-x64-linux.tar.gz
ENV PATH="/opt/ncbi-blast-2.17.0+/bin:${PATH}"

# ── Reference proteome FASTA (Ensembl GRCh38) ────────────────────────────────
RUN mkdir -p /opt/ref_proteome \
    && wget -q https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz \
        -O /opt/ref_proteome/Homo_sapiens.GRCh38.pep.all.fa.gz \
    && gunzip /opt/ref_proteome/Homo_sapiens.GRCh38.pep.all.fa.gz

# set permissions
RUN chown -R $MAMBA_USER:$MAMBA_USER /opt/PRIME /opt/MixMHCpred /opt/ncbi-blast-2.17.0+ /opt/ref_proteome
# change user
USER $MAMBA_USER

# set workdir
WORKDIR /home/ec2-user

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

