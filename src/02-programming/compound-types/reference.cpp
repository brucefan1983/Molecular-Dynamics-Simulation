#include <iostream> // std::cout, std::endl;

int main()
{
  int a = 1;
  int& ra = a;
  std::cout << "a = " << a << std::endl;   // 1
  std::cout << "ra = " << ra << std::endl; // 1
  ra = 2;
  std::cout << "a = " << a << std::endl;   // 2
  std::cout << "ra = " << ra << std::endl; // 2
  return 0;
}
