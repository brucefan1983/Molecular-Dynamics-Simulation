#include <iostream> // std::cout, std::endl;

int main()
{
  int n = 6;
  int f = 1;
  for (int i = 1; i <= n; ++i) {
    f *= i;
    std::cout << i << "! = " << f << std::endl;
  }

  return 0;
}