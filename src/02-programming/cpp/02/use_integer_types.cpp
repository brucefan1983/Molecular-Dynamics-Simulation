#include <iostream>

int main()
{
  bool fact = false;
  std::cout << fact << std::endl; // 0
  fact = true;
  std::cout << fact << std::endl; // 1

  char a_character = 'a';
  std::cout << a_character << std::endl; // a
  a_character = '8';
  std::cout << a_character << std::endl; // 8

  int an_integer = -3;
  std::cout << an_integer << std::endl; // -3
  int another_integer = an_integer;
  std::cout << another_integer << std::endl; // -3

  return 0;
}