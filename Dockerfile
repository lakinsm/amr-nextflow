FROM ubuntu:15.10

MAINTAINER Steven Lakin <Steven.Lakin@colostate.edu>

RUN git clone https://github.com/cdeanj/coverage_sampler.git && \
	cd coverage_sampler && \
	make && \
	cp csa ${HOME}/bin

