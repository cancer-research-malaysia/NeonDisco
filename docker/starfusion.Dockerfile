# bring in the micromamba image so we can copy files from it
FROM mambaorg/micromamba as micromamba

# This is the image we are going add micromaba to: #star-fusion version 1.15.0, fusioninspector version 2.10.0
FROM trinityctat/starfusion 

LABEL maintainer="Suffian Azizan"
LABEL version="3.0"
LABEL description="container image of STAR Fusion from TrinityCTAT"

# change to root user
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

# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed less nano \
&& rm -rf /var/lib/apt/lists/*

## copy files
COPY starfusion/src/fusion_lib.Mar2021.dat.gz /tmp/fusion_lib.Mar2021.dat.gz
COPY starfusion/src/AnnotFilterRule.pm /tmp/AnnotFilterRule.pm

# delete source codes that require hotfixes
RUN rm -f /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl && rm -f /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/util/mask_confounding_features_from_genome.pl
#replace the default prep_genome_lib.pl script with a hotfix for PAR region by github/zhuchcn (SHA: dc5c779a21f23f5c07e67850ba6845da7b058234)
COPY --chown=root:root starfusion/src/prep_genome_lib-par-hotfix.pl /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl
COPY --chown=root:root starfusion/src/mask_confounding_features_from_genome-par-hotfix.pl /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/util/mask_confounding_features_from_genome.pl
RUN chmod +x /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl && \
    chmod +x /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/util/mask_confounding_features_from_genome.pl

# copy edited FusionInspector script
COPY --chown=root:root starfusion/src/FusionInspector /tmp/FusionInspector
# delete the original FusionInspector script and replace it with the edited one
RUN rm /usr/local/src/STAR-Fusion/FusionInspector/FusionInspector && mv /tmp/FusionInspector /usr/local/src/STAR-Fusion/FusionInspector/
# make the FusionInspector script executable
RUN chmod +x /usr/local/src/STAR-Fusion/FusionInspector/FusionInspector

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
COPY --chown=$MAMBA_USER:$MAMBA_USER starfusion/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# set workdir
WORKDIR /home/ec2-user

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
