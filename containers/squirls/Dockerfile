FROM ubuntu:22.04

# Configuring TZ to avoid halting git installation
ENV TZ=Europe/London
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Installing necessary utilities and dependencies
RUN  apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install -y default-jdk \
	&& apt-get install -y git-all

RUN git clone https://github.com/TheJacksonLaboratory/Squirls \
	&& cd Squirls \
	&& ./mvnw package