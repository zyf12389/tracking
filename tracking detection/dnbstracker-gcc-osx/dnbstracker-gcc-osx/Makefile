CC = clang++
CPPFLAGS = -c -Wall -O3 -DNDEBUG -I/stow/opencv-2.4.10/include -I/stow/opencv-2.4.10/include/opencv -I/usr/include/libxml2
LDFLAGS = -O3 $(shell pkg-config --libs opencv) /usr/lib/libxml2.dylib
SOURCES = DnbsClassAlpha.cpp DnbsTracker.cpp    NaiveNbsClass.cpp  utilities.cpp DnbsClassBeta.cpp  HaarBase.cpp       commonlib.cpp DnbsClassGamma.cpp LyonLib.cpp        main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = dnbstrack

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o: 
	$(CC) $(CPPFLAGS) $< -o $@

clean:
	rm -rf *.o
	rm -rf $(EXECUTABLE)
