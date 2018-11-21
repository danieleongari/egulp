include gcc.arch

all: egulp
egulp: egulp.o geometry.o scf.o utilsa.o genpot.o interact.o param.o vdw.o grid.o cell.o econfig.o pcharge.o eseqparam.o psinv.o elements.o 
	${CC} ${CFLAGS}  egulp.o geometry.o scf.o utilsa.o genpot.o interact.o param.o vdw.o grid.o cell.o econfig.o pcharge.o eseqparam.o psinv.o elements.o `pkg-config --cflags --libs gsl` -o egulp
egulpomp: egulp.o geometry.o scf.o utilsa.o genpotomp.o interact.o param.o vdw.o grid.o cell.o econfig.o pcharge.o eseqparam.o psinv.o elements.o 
	${CC} ${CFLAGS}  egulp.o geometry.o scf.o utilsa.o genpotomp.o interact.o param.o vdw.o grid.o cell.o econfig.o pcharge.o eseqparam.o psinv.o elements.o `pkg-config --cflags --libs gsl` -openmp -lpthread -o egulpomp
egulp.o: egulp.c geometry.h scf.h utilsa.h genpot.h param.h grid.h econfig.h pcharge.h eseqparam.h elements.h  
	${CC} ${CFLAGS}  -c egulp.c -o egulp.o
geometry.o: geometry.c geometry.h utilsa.h vdw.h cell.h
	${CC} ${CFLAGS} -c geometry.c -o geometry.o
vdw.o: vdw.c
	${CC} ${CFLAGS} -c vdw.c -o vdw.o	
grid.o: grid.c grid.h utilsa.h vdw.h 
	${CC} ${CFLAGS} -c grid.c -o grid.o
cell.o: cell.c 
	${CC} ${CFLAGS} -c cell.c -o cell.o	
scf.o: scf.c scf.h utilsa.h
	${CC} ${CFLAGS} -c scf.c -o scf.o	
utilsa.o: utilsa.c 
	${CC} ${CFLAGS}  -c utilsa.c -o utilsa.o
genpot.o: genpot.c utilsa.h geometry.h interact.h cell.h elements.h
	${CC} ${CFLAGS} -c genpot.c -o genpot.o
interact.o: interact.c
	${CC} ${CFLAGS} -c interact.c -o interact.o
param.o: param.c elements.h
	${CC} ${CFLAGS} -c param.c -o param.o
econfig.o: econfig.c econfig.h
	${CC} ${CFLAGS} -c econfig.c -o econfig.o
pcharge.o: pcharge.c pcharge.h utilsa.h
	${CC} ${CFLAGS} -c pcharge.c -o pcharge.o
genpotomp.o: genpotomp.c utilsa.h geometry.h interact.h cell.h elements.h
	${CC} ${CFLAGS}  -openmp -c genpotomp.c -o genpotomp.o
eseqparam.o: eseqparam.c eseqparam.h utilsa.h  psinv.h elements.h
	${CC} ${CFLAGS} -c eseqparam.c -o eseqparam.o
psinv.o:
	${CC} ${CFLAGS} -c psinv.c -o psinv.o
elements.o:
	${CC} ${CFLAGS} -c elements.c -o elements.o	
clean:
	rm -f egulp.o geometry.o scf.o utilsa.o genpot.o interact.o param.o vdw.o cell.o grid.o econfig.o pcharge.o genpotomp.o eseqparam.o psinv.o elements.o
