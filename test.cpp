#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main (int argc, char* argv[])
{
  if (argc != 3)
  {
    printf("fail\n");
  }
  else 
  {
    int num1 = atoi(argv[1]);
    int num2 = atoi(argv[2]);
    std::cout << num2/num1 << std::endl;
  }
}