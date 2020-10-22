#include <iostream>

int main()
{
  char oper;
  float x = 1.0;
  float y = 2.0;
  std::cout << "Enter an operator (+, -, *, /): ";
  std::cin >> oper;

  switch (oper) {
    case '+':
      std::cout << x << " + " << y << " = " << x + y;
      break;
    case '-':
      std::cout << x << " - " << y << " = " << x - y;
      break;
    case '*':
      std::cout << x << " * " << y << " = " << x * y;
      break;
    case '/':
      std::cout << x << " / " << y << " = " << x / y;
      break;
    default:
      std::cout << "Error: Invalid operator.";
      break;
  }

  return 0;
}