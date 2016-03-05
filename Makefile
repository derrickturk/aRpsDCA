VER=1.0.3

aRpsDCA_$(VER).tar.gz:
	R CMD build .

install: aRpsDCA_$(VER).tar.gz uninstall
	R CMD INSTALL aRpsDCA_$(VER).tar.gz

clean:
	-rm aRpsDCA_$(VER).tar.gz

uninstall:
	-R CMD REMOVE aRpsDCA

check: aRpsDCA_$(VER).tar.gz
	R CMD check aRpsDCA_$(VER).tar.gz
