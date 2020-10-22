#include <iostream> // std::cout, std::endl;

int main()
{
  const int n = 6;
  int f = 1;
  int i = 1;
  while (true) {
    f *= i;
    if (f > 1000)
      break;
    std::cout << i << "! = " << f << std::endl;
    ++i;
  }

  return 0;
}
