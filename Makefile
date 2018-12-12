all: strassen.cpp
	g++ -o strassen strassen.cpp

strassen: strassen.cpp
	g++ -o strassen strassen.cpp

clean:
	$(RM) strassen