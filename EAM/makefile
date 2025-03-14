objectes = AwadLRO.o Belashchenko.o Zhou.o
compilador = gfortran
OPT = -O3

eam.file.x: $(objectes)
	$(compilador) main.f90 -o eam.file.x $(OPT) $(objectes)
arrels.x: Zhou.o AwadLRO.o Belashchenko.o
	$(compilador) ARRELS.f90 -o arrels.x $(OPT) Zhou.o AwadLRO.o Belashchenko.o

AwadLRO.o: AwadLRO.f90
	$(compilador) -c $(OPT) AwadLRO.f90
Belashchenko.o: Belashchenko.f90
	$(compilador) -c $(OPT) Belashchenko.f90
Zhou.o: Zhou.f90
	$(compilador) -c $(OPT) Zhou.f90


.PHONY: clean move
clean:
	rm -f *.o
	rm -f *mod
	rm -f *eam.alloy
move:
	mv *eam.alloy EAM.FILES