# See ../triqs/packaging for other options
FROM flatironinstitute/triqs:unstable-ubuntu-clang
ARG APPNAME=keldy

COPY requirements.txt /src/$APPNAME/requirements.txt
RUN pip install -r /src/$APPNAME/requirements.txt

ADD https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz /tmp/
RUN cd /tmp && tar xf gsl-2.6.tar.gz && cd gsl-2.6 && \
      ./configure --prefix=/opt/gsl-2.6 && make -j${NTHREAD} && make install && \
      cd /tmp && rm -rf /tmp/gsl-2.6*
ENV GSL_ROOT=/opt/gsl-2.6  

COPY --chown=build . $SRC/$APPNAME
WORKDIR $BUILD/$APPNAME
RUN chown build .
USER build
ARG BUILD_DOC=0
RUN cmake $SRC/$APPNAME -DTRIQS_ROOT=${INSTALL} -DBuild_Documentation=${BUILD_DOC} && make -j2 && make test CTEST_OUTPUT_ON_FAILURE=1
USER root
RUN make install
