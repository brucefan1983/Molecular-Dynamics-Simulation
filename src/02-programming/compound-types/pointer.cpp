#include <iostream> // std::cout, std::endl;

int main()
{
  int a = 1;
  int* pa = &a;
  std::cout << "a = " << a << std::endl;     // value of a
  std::cout << "&a = " << &a << std::endl;   // address of a
  std::cout << "pa = " << pa << std::endl;   // address of a
  std::cout << "*pa = " << *pa << std::endl; // value of a
  return 0;
}
