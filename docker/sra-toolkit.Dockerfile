FROM ncbi/sra-tools

LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="minimal container image for SRA Toolkit for CRMY"

# Switch to root for installation
USER root

# Install MatchHostFsOwner. Using version 1.0.1
# See https://github.com/FooBarWidget/matchhostfsowner/releases
ADD https://github.com/FooBarWidget/matchhostfsowner/releases/download/v1.0.1/matchhostfsowner-1.0.1-x86_64-linux.gz /sbin/matchhostfsowner.gz
RUN gunzip /sbin/matchhostfsowner.gz && \
  chown root: /sbin/matchhostfsowner && \
  chmod +x /sbin/matchhostfsowner

# Create group and user with Alpine/BusyBox syntax
RUN addgroup -g 9999 app && \
  adduser -D -u 9999 -G app -s /bin/sh -h /home/app app

# Set workdir
WORKDIR /home/app

# Create the entrypoint directory if it doesn't exist
RUN mkdir -p /usr/local/bin

# Create the entrypoint script
RUN echo '#!/bin/sh' > /usr/local/bin/_entrypoint.sh && \
    echo 'exec "$@"' >> /usr/local/bin/_entrypoint.sh && \
    chmod +x /usr/local/bin/_entrypoint.sh

# Verify the entrypoint script exists and is executable
RUN ls -la /usr/local/bin/_entrypoint.sh && \
    cat /usr/local/bin/_entrypoint.sh

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]