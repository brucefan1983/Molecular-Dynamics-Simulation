## 程序 simpleMD.cpp


该 C++ 程序可以在 Linux 和 Windows 操作系统使用。我们推荐使用 GCC 工具。在命令行可以用如下方式编译本章的程序：
```shell
    $ g++ -O3 simpleMD.cpp -o simpleMD
```
其中，`-O3` 选项表示优化等级。编译完成后，将生成名为 `simpleMD` 的可执行文件（在 Windows 中为 `simpleMD.exe`)。

然后，就可以在命令行使用该程序：
```shell
    $ simpleMD numCells numSteps temperature timeStep
```
如果使用时忘了给命令行参数，程序会提示正确的用法。

具体地，笔者用如下命令运行程序：
```shell
    $ simpleMD 4 20000 60 5
```
也就是说，考虑原子数为 $4^3 \times 4 = 256$ 的固态氩体系，运行 $2\times 10^4$ 步，初始温度为 60 K，积分步长为 5 fs。该模拟在笔者的计算机中运行的时间约为 10 秒钟。


