/*
    Copyright 2017 Zheyong Fan
    This file is part of XiaoMD.
    XiaoMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    XiaoMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with XiaoMD.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "xiaomd.h"
#include <ctime>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{

  const clock_t tStart = clock();

  XiaoMD xiaomd;

  const clock_t tStop = clock();
  const float tElapsed = float(tStop - tStart) / CLOCKS_PER_SEC;
  std::cout << "Time used = " << tElapsed << " s" << std::endl;

  return 0;
}
