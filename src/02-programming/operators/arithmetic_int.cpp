#include <iostream> // std::cout, std::endl;

int main()
{
  int x = 10;
  int y = -3;
  int z = -10;

  std::cout << x + y << std::endl; // 7
  std::cout << x - y << std::endl; // 13
  std::cout << x * y << std::endl; // -30

  std::cout << x / y << std::endl;                 // -3
  std::cout << z / y << std::endl;                 // 3
  std::cout << x % y << std::endl;                 // 1
  std::cout << z % y << std::endl;                 // -1
  std::cout << (x / y) * y + (x % y) << std::endl; // 10
  std::cout << (z / y) * y + (z % y) << std::endl; // -10

  int a = ++x;
  std::cout << a << std::endl; // 11
  std::cout << x << std::endl; // 11

  a = x++;
  std::cout << a << std::endl; // 11
  std::cout << x << std::endl; // 12

  --x;
  std::cout << x << std::endl; // 11

  x--;
  std::cout << x << std::endl; // 10

  return 0;
}