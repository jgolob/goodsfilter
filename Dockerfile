# golob/goodsfilter:0.1.0
#

FROM      ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
python3-dev \
python3-pip \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install pip --upgrade
RUN pip3 install \
awscli>=1.15.14 \
boto3>=1.7.14
ADD goodsfilter/goodsfilter.py /usr/local/bin/goodsfilter
RUN chmod +x /usr/local/bin/goodsfilter
