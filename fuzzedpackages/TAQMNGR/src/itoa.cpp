#include <taq.h>


void itoa(int a, char *ChOut)
{
	int x, b = a;
	char i;
    while(a)
    {
        x = a%10;
        a /= 10;
        ChOut++;
    }
    *ChOut = '\0';
    while(b)
    {
		x = b%10;
		b /= 10;
		i = '0';
		i = i+x;
		ChOut--;
		*ChOut = i;
	}
}
