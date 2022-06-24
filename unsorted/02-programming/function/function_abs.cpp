#include <iostream> // std::cout, std::endl;

float my_abs1(float x) // can modify the parameter x
{
  if (x < 0.0f)
    x = -x;
  return x;
}

float my_abs2(const float x) // cannot modify the parameter x
{
  float y = x;
  if (x < 0.0f)
    y = -x;
  return y;
}

float my_abs3(const float x) // can return multiple times
{
  if (x < 0.0f)
    return -x;
  else
    return x;
}

float my_abs4(const float x); // declaration
float my_abs4(const float x); // can declare more than once
float my_abs4(const float y); // parameter name can be changed
float my_abs4(const float);   // parameter name can be omittd

int main()
{
  float x = -1.0f;
  std::cout << "my_abs1(x) = " << my_abs1(x) << std::endl;
  std::cout << "my_abs2(x) = " << my_abs2(x) << std::endl;
  std::cout << "my_abs3(x) = " << my_abs3(x) << std::endl;
  std::cout << "my_abs4(x) = " << my_abs4(x) << std::endl;
  return 0;
}

float my_abs4(const float x)
{
  if (x < 0.0f)
    return -x;
  else
    return x;
}
