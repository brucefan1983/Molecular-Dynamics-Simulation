#include <iostream>

int main()
{
  char oper;
  float x = 1.0;
  float y = 2.0;
  std::cout << "Enter an operator (+, -, *, /): ";
  std::cin >> oper;

  if (oper == '+') {
    std::cout << x << " + " << y << " = " << x + y;
  } else if (oper == '-') {
    std::cout << x << " - " << y << " = " << x - y;
  } else if (oper == '*') {
    std::cout << x << " * " << y << " = " << x * y;
  } else if (oper == '/') {
    std::cout << x << " / " << y << " = " << x / y;
  } else {
    std::cout << "Error: Invalid operator.";
  }

  return 0;
}