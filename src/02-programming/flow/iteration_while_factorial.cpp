#include <iostream> // std::cout, std::endl;

int main()
{
  const int n = 6;
  int f = 1;
  int i = 1;
  while (i <= 6) {
    f *= i;
    std::cout << i << "! = " << f << std::endl;
    ++i;
  }

  return 0;
}