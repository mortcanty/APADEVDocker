# Acquisition Path Analysis

This repo contains the source files for the Docker image mort/apadevdocker, which performs acquisition path analysis on a hypothetical nuclear fuel cycle.
###Installation

 1. <a href="https://docs.docker.com/installation/ubuntulinux/">install Docker</a>
 2. In a terminal, run the command<br />
      __sudo docker run -d -p 433:8888 --name=apa  mort/apadevdocker__<br />
 3. Point your browser to<br /> 
    __localhost:433__
 4. Open a new notebook with __New/Python 2__.
 5. Test with __%run pm__