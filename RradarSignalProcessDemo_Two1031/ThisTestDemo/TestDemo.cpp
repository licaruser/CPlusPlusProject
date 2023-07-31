#include <malloc.h>
#include <stdio.h>


int* remove_element(int* array, int sizeOfArray, int indexToRemove)
{
	int* temp = (int*)malloc((sizeOfArray - 1) * sizeof(int)); // allocate an array with a size 1 less than the current one if (indexToRemove != 0) 
	memcpy(temp, array, indexToRemove * sizeof(int)); // copy everything BEFORE the index if (indexToRemove != (sizeOfArray - 1)) memcpy(temp+indexToRemove, array+indexToRemove+1, 
	(sizeOfArray - indexToRemove - 1) * sizeof(int)); // copy everything AFTER the index 
	free(array);
	return temp;
}
int main()
{
	int howMany = 20;
	int* test = (int*)malloc(howMany * sizeof(int*));
	for (int i = 0; i < howMany; ++i)
		(test[i]) = i;
	printf("%dn", test[16]);
	remove_element(test, howMany, 16);
	--howMany;
	printf("%dn", test[16]); return 0;
}