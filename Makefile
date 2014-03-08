VER=0.1.0

aRpsDCA_$(VER).tar.gz:
	R CMD BUILD .

install: aRpsDCA_$(VER).tar.gz uninstall
	R CMD INSTALL aRpsDCA_$(VER).tar.gz

clean:
	-rm aRpsDCA_$(VER).tar.gz

uninstall:
	-R CMD REMOVE aRpsDCA

check: aRpsDCA_$(VER).tar.gz
	R CMD CHECK aRpsDCA_$(VER).tar.gz
