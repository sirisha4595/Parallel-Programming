INC=/lib/modules/$(shell uname -r)/build/arch/x86/include

all: sample rArB

sample: sample.cpp
	mpicc -Wall -I$(INC)/generated/uapi -I$(INC)/uapi sample.cpp -o sample

rArB: rArB.cpp
	mpicc -Wall -I$(INC)/generated/uapi -I$(INC)/uapi rArB.cpp -o rArB

clean:
	make -C /lib/modules/$(shell uname -r)/build M=$(PWD) clean
	rm -f sample

