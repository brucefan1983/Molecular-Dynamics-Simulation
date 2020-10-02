#include <iostream> // std::cout, std::endl;

int main()
{
  double x = 1.2;
  double y = 3.0;

  std::cout << x + y << std::endl;         // 4.2
  std::cout << x - y << std::endl;         // -1.8
  std::cout << x * y << std::endl;         // 3.6
  std::cout << x / y << std::endl;         // 0.4
  std::cout << x + y * 2.0 << std::endl;   // 7.2
  std::cout << (x + y) * 2.0 << std::endl; // 8.4
  std::cout << x / y * 2.0 << std::endl;   // 0.8
  std::cout << x / (y * 2.0) << std::endl; // 0.2

  return 0;
}