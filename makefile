all: main.exe

main.exe: runscript.cpp
	g++ -o runscript.exe runscript.cpp
	./runscript.exe

runscript.exe: template.cpp
	g++ -o template.exe template.cpp
	./template.exe