FROM griffithlab/pvactools:latest
LABEL maintainer="Suffian Azizan"
LABEL version="3.0"
LABEL description="PVACtools with MatchHostFsOwner for rootless Docker compatibility"

# Switch to root to install MatchHostFsOwner
USER root

# Install MatchHostFsOwner. Using version 1.0.1
# See https://github.com/FooBarWidget/matchhostfsowner/releases
ADD https://github.com/FooBarWidget/matchhostfsowner/releases/download/v1.0.1/matchhostfsowner-1.0.1-x86_64-linux.gz /sbin/matchhostfsowner.gz
RUN gunzip /sbin/matchhostfsowner.gz && \
  chown root: /sbin/matchhostfsowner && \
  chmod +x /sbin/matchhostfsowner

# Create the app user/group that MatchHostFsOwner will use
RUN addgroup --gid 9999 app && \
  adduser --uid 9999 --gid 9999 --disabled-password --gecos App app

# Set working directory
WORKDIR /home/app

# Since the original image just uses CMD ["/bin/bash"], we need a simple entrypoint
ENTRYPOINT ["/sbin/matchhostfsowner"]
CMD ["/bin/bash"]

