FROM images.sbgenomics.com/aberno/rhapsody:1.8

RUN echo "source ~/.bash_profile" >> /root/.bashrc
COPY .bash_profile /root

RUN apt-get install -y -qq --no-install-recommends vim