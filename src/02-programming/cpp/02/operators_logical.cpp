#include <iostream> // std::cout, std::fixed, std::endl;

int main()
{
  bool relation = true && true;       // true
  std::cout << relation << std::endl; // 1
  relation = true && false;           // false
  std::cout << relation << std::endl; // 0
  relation = false && true;           // false
  std::cout << relation << std::endl; // 0
  relation = false && false;          // false
  std::cout << relation << std::endl; // 0

  relation = true || true;            // true
  std::cout << relation << std::endl; // 1
  relation = true || false;           // true
  std::cout << relation << std::endl; // 1
  relation = false || true;           // true
  std::cout << relation << std::endl; // 1
  relation = false || false;          // false
  std::cout << relation << std::endl; // 0

  relation = !true;                   // false
  std::cout << relation << std::endl; // 0
  relation = !false;                  // true
  std::cout << relation << std::endl; // 1

  return 0;
}