VER=0.1.0

install: aRpsDCA_$(VER).tar.gz clean
	R CMD INSTALL aRpsDCA_$(VER).tar.gz

aRpsDCA_$(VER).tar.gz:
	R CMD BUILD .

clean:
	R CMD REMOVE aRpsDCA
