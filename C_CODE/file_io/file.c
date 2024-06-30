// C Program to demonstrate fscanf
#include <stdio.h>

// Driver Code
int main()
{
    FILE* ptr = fopen("info.txt", "r");
    if (ptr == NULL) {
        printf("no such file.");
        return 0;
    }

    /* Assuming that abc.txt has content in below
       format
       NAME    AGE   CITY
       abc     12    hyderabad
       bef     25    delhi
       cce     65    bangalore */
    char name[100];
    int age;
    char city[100];
    while (fscanf(ptr, "%s %d %s ", name, &age, city) != EOF)
        printf("%s %d %s\n", name, age, city);

    return 0;
}