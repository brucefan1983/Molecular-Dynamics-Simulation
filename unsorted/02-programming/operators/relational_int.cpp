#include <iostream> // std::cout, std::endl;

int main()
{
  int a = 10;
  int b = -3;
  int c = 20;
  int d = a; // d = 10

  bool relation = (a == b); // false
  std::cout << relation << std::endl;
  relation = (a == d); // true
  std::cout << relation << std::endl;

  relation = (a != b); // true
  std::cout << relation << std::endl;
  relation = (a != d); // false
  std::cout << relation << std::endl;

  relation = (a > b); // true
  std::cout << relation << std::endl;
  relation = (a > c); // false
  std::cout << relation << std::endl;

  relation = (a < b); // false
  std::cout << relation << std::endl;
  relation = (a < c); // true
  std::cout << relation << std::endl;

  relation = (a >= c); // false
  std::cout << relation << std::endl;
  relation = (a >= d); // true
  std::cout << relation << std::endl;

  relation = (a <= c); // true
  std::cout << relation << std::endl;
  relation = (a <= d); // true
  std::cout << relation << std::endl;

  return 0;
}