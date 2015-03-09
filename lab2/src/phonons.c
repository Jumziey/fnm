#include <stdio.h>
#include <stdlib.h>
//#include "frequencies.c"
int
main(int argc, char** argv)
{
	int i;
	for(i=0;i<3;i++)
		printf("%d %d %d\n", i%3, (i+1)%3, (i+2)%3);
	return 0;
}
